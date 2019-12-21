#include "minacalc.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <memory>
#include <numeric>
#include <iostream>

using std::cout;
using std::endl;

using std::vector;
using std::min;
using std::max;
using std::sqrt;
using std::pow;

template<typename T>
T CalcClamp(T x, T l, T h) {
return x > h ? h : (x < l ? l : x);
}

template <typename F>
float approximate(float value, float resolution, int num_iters, F is_too_low, bool limit_at_100 = false) {
    for (int i = 0; i < num_iters; i++) {
        while (is_too_low(value)) {
            if (limit_at_100 && value > 100.f) return value;
            value += resolution;
        }
        value -= resolution;
        resolution /= 2.f;
    }
    
    return value + 2.f * resolution;
}

inline float mean(const vector<float>& v) {
    return std::accumulate(begin(v), end(v), 0.f) / v.size();
}

// Coefficient of variation
inline float cv(const vector<float> &input) {
    float sd = 0.f;
    float average = mean(input);
    for (float i : input)
        sd += (i - average)*(i - average);

    return sqrt(sd / input.size()) / average;
}

inline float downscale_low_accuracy_scores(float f, float sg) {
    return sg >= 0.93f ? f : min(max(f - sqrt(0.93f - sg), 0.f), 100.f);
}

// Moving average with n=3. The `neutral` value is used for the
// "out-of-bounds values" required for the moving averages on the start
// and end.
inline void Smooth(vector<float>& input, float neutral) {
    float f1;
    float f2 = neutral;
    float f3 = neutral;

    for (float & i : input) {
        f1 = f2;
        f2 = f3;
        f3 = i;
        i = (f1 + f2 + f3) / 3;
    }
}

// Like `Smooth()`, but with n=2 and neutral value zero.
inline void DifficultyMSSmooth(vector<float>& input) {
    float f1;
    float f2 = 0.f;

    for (float & i : input) {
        f1 = f2;
        f2 = i;
        i = (f1 + f2) / 2.f;
    }
}

inline float AggregateScores(const vector<float>& skillsets, float rating, float resolution) {
    auto check_if_too_low = [skillsets](float rating) {
        float sum = 0.0f;
        for (float i : skillsets) {
            sum += 2.f / std::erfc(0.5f * (i - rating)) - 1.f;
        }
        return 3 < sum;
    };
    return approximate(rating, resolution, 11, check_if_too_low);
}

// Converts a row byte into the number of taps present in the row
// e.g. 0010 -> 1 or 1011 -> 3
unsigned int column_count(unsigned int note) {
    return note % 2 + note / 2 % 2 + note / 4 % 2 + note / 8 % 2;
}

// Proportion of how many chords of size `chord_size` exist
float chord_proportion(const vector<NoteInfo>& NoteInfo, const int chord_size) {
    unsigned int taps = 0;
    unsigned int chords = 0;

    for (auto row : NoteInfo) {
        unsigned int notes = column_count(row.notes);
        taps += notes;
        if (notes == chord_size)
            chords += notes;
    }

    return static_cast<float>(chords) / static_cast<float>(taps);
}

vector<float> skillset_vector(const DifficultyRating& difficulty) {
    return vector<float> {difficulty.overall,
                          difficulty.stream,
                          difficulty.jumpstream,
                          difficulty.handstream,
                          difficulty.stamina,
                          difficulty.jack,
                          difficulty.chordjack,
                          difficulty.technical
    };
}

float highest_difficulty(const DifficultyRating& difficulty) {
    auto v = {difficulty.stream,difficulty.jumpstream,difficulty.handstream,difficulty.stamina,difficulty.jack,
              difficulty.chordjack,difficulty.technical};
    return *std::max_element(v.begin(), v.end());
}

DifficultyRating Calc::CalcMain(const vector<NoteInfo>& NoteInfo, float music_rate, float score_goal) {
    float grindscaler = CalcClamp(0.93f + (0.07f * (NoteInfo.back().rowTime - 30.f) / 30.f), 0.93f, 1.f)
            * CalcClamp(0.873f + (0.13f * (NoteInfo.back().rowTime - 15.f) / 15.f), 0.87f, 1.f);

    float shortstamdownscaler = CalcClamp(0.9f + (0.1f * (NoteInfo.back().rowTime - 150.f) / 150.f), 0.9f, 1.f);

    float jprop = chord_proportion(NoteInfo, 2);
    float nojumpsdownscaler = CalcClamp(0.8f + (0.2f * (jprop + 0.5f)), 0.8f, 1.f);
    float manyjumpsdownscaler = CalcClamp(1.43f - jprop, 0.85f, 1.f);

    float hprop = chord_proportion(NoteInfo, 3);
    float nohandsdownscaler = CalcClamp(0.8f + (0.2f * (hprop + 0.75f)), 0.8f, 1.f);
    float allhandsdownscaler = CalcClamp(1.23f - hprop, 0.85f, 1.f);

    float qprop = chord_proportion(NoteInfo, 4);
    float lotquaddownscaler = CalcClamp(1.13f - qprop, 0.85f, 1.f);

    float jumpthrill = CalcClamp(1.625f - jprop - hprop, 0.85f, 1.f);

    InitializeHands(NoteInfo, music_rate);
    TotalMaxPoints();
    float stream = Chisel(0.1f, 10.24f, score_goal, CHISEL_NPS);
    float js = Chisel(0.1f, 10.24f, score_goal, CHISEL_NPS | CHISEL_JS);
    float hs = Chisel(0.1f, 10.24f, score_goal, CHISEL_NPS | CHISEL_HS);
    float tech = Chisel(0.1f, 10.24f, score_goal, 0);
    float jack = Chisel(0.1f, 10.24f, score_goal, CHISEL_NPS | CHISEL_JACK);

    float techbase = max(stream, jack);
    tech *= CalcClamp(tech / techbase, 0.85f, 1.f);

    float stam;
    if (stream > tech || js > tech || hs > tech)
        if (stream > js && stream > hs)
            stam = Chisel(stream - 0.1f, 2.56f, score_goal, CHISEL_STAM | CHISEL_NPS);
        else if (js > hs)
            stam = Chisel(js - 0.1f, 2.56f, score_goal, CHISEL_STAM | CHISEL_NPS | CHISEL_JS);
        else
            stam = Chisel(hs - 0.1f, 2.56f, score_goal, CHISEL_STAM | CHISEL_NPS | CHISEL_HS);
    else
        stam = Chisel(tech - 0.1f, 2.56f, score_goal, CHISEL_STAM);

    js *= 0.95f;
    hs *= 0.95f;
    stam *= 0.9f;

    float chordjack = jack * 0.75f;
    tech *= 0.95f;

    DifficultyRating difficulty = DifficultyRating {0.0,
                                                    downscale_low_accuracy_scores(stream, score_goal),
                                                    downscale_low_accuracy_scores(js, score_goal),
                                                    downscale_low_accuracy_scores(hs, score_goal),
                                                    downscale_low_accuracy_scores(stam, score_goal),
                                                    downscale_low_accuracy_scores(jack, score_goal),
                                                    downscale_low_accuracy_scores(chordjack, score_goal),
                                                    downscale_low_accuracy_scores(tech, score_goal)
    };

    chordjack = difficulty.handstream;

    difficulty.stream *= allhandsdownscaler * manyjumpsdownscaler * lotquaddownscaler;
    difficulty.jumpstream *= nojumpsdownscaler * allhandsdownscaler * lotquaddownscaler;
    difficulty.handstream *= nohandsdownscaler * allhandsdownscaler * 1.015f * manyjumpsdownscaler * lotquaddownscaler;
    difficulty.stamina = CalcClamp(difficulty.stamina * shortstamdownscaler * 0.985f * lotquaddownscaler, 1.f,
                                   max(max(difficulty.stream, difficulty.jack), max(difficulty.jumpstream, difficulty.handstream)) * 1.1f);
    difficulty.technical *= allhandsdownscaler * manyjumpsdownscaler * lotquaddownscaler * 1.01f;

    chordjack *= CalcClamp(qprop + hprop + jprop + 0.2f, 0.5f, 1.f) * 1.025f;

    bool downscale_chordjack_at_end = false;
    if (chordjack > difficulty.jack)
        difficulty.chordjack = chordjack;
    else
        downscale_chordjack_at_end = true;

    fingerbias /= static_cast<float>(2 * nervIntervals.size());
    float finger_bias_scaling = CalcClamp(3.55f - fingerbias, 0.85f, 1.f);
    difficulty.technical *= finger_bias_scaling;

    if (finger_bias_scaling <= 0.95f) {
        difficulty.jack *= 1.f + (1.f - sqrt(finger_bias_scaling));
    }
    
    // If HS or JS are more prominent than stream, downscale stream a
    // little to prevent too much stream rating as a side effect from
    // JS/HS.
    // Stream is nerfed by `sqrt(hs - stream)` or `sqrt(js - stream)`
    float max_js_hs = max(difficulty.handstream, difficulty.jumpstream);
    if (difficulty.stream < max_js_hs)
        difficulty.stream -= sqrt(max_js_hs - difficulty.stream);

    vector<float> temp_vec = skillset_vector(difficulty);
    float overall = AggregateScores(temp_vec, 0.f, 10.24f);
    difficulty.overall = downscale_low_accuracy_scores(overall, score_goal);

    temp_vec = skillset_vector(difficulty);
    float aDvg = mean(temp_vec) * 1.2f;
    difficulty.overall = downscale_low_accuracy_scores(min(difficulty.overall, aDvg) * grindscaler, score_goal);
    difficulty.stream = downscale_low_accuracy_scores(min(difficulty.stream, aDvg * 1.0416f) * grindscaler, score_goal);
    difficulty.jumpstream = downscale_low_accuracy_scores(min(difficulty.jumpstream, aDvg * 1.0416f) * grindscaler, score_goal) * jumpthrill;
    difficulty.handstream = downscale_low_accuracy_scores(min(difficulty.handstream, aDvg) * grindscaler, score_goal) * jumpthrill;
    difficulty.stamina = downscale_low_accuracy_scores(min(difficulty.stamina, aDvg) * grindscaler, score_goal) * sqrt(jumpthrill) * 0.996f;
    difficulty.jack = downscale_low_accuracy_scores(min(difficulty.jack, aDvg) * grindscaler, score_goal);
    difficulty.chordjack = downscale_low_accuracy_scores(min(difficulty.chordjack, aDvg) * grindscaler, score_goal);
    difficulty.technical = downscale_low_accuracy_scores(min(difficulty.technical, aDvg * 1.0416f) * grindscaler, score_goal) * sqrt(jumpthrill);

    float highest = max(difficulty.overall, highest_difficulty(difficulty));

    vector<float> temp = skillset_vector(difficulty);
    difficulty.overall = AggregateScores(temp, 0.f, 10.24f);

    if (downscale_chordjack_at_end) {
        difficulty.chordjack *= 0.9f;
    }

    // Calculate and check minimum required percentage. This percentage
    // is dependant on MSD value. It's a linear function, clamped
    // between 50% and 90%. It starts at `0 MSD -> 50%` and ends at
    // `40 MSD -> 90%`
    float minimum_required_percentage = CalcClamp(0.5f + (highest / 100.f), 0.f, 0.9f);
    if (score_goal < minimum_required_percentage) {
        difficulty = DifficultyRating {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
    }

    // If technical is supposedly the highest skillset, but JS or HS are
    // near to it, technical might be falsely rated too high. In that
    // case downscale
    if (highest == difficulty.technical) {
        auto hs = difficulty.handstream;
        auto js = difficulty.jumpstream;
        
        // If technical within 4.5 points of HS or JS, downscale it.
        difficulty.technical -= CalcClamp(4.5f - (difficulty.technical - hs), 0.f, 4.5f);
        difficulty.technical -= CalcClamp(4.5f - (difficulty.technical - js), 0.f, 4.5f);
    }

    difficulty.jack *= 0.925f;
    difficulty.technical *= 1.025f;
    difficulty.overall = highest_difficulty(difficulty);

    return difficulty;
}

int jack_loss_n = 0;
// ugly jack stuff
// Parameter j = JackSeq
// Parameter x = potential player rating
float Calc::JackLoss(const vector<float>& j, float x) {
    float output = 0.f;
    float ceiling = 1.f;
    float mod = 1.f;
    float base_ceiling = 1.15f;
    float fscale = 1750.f;
    float prop = 0.75f;
    float mag = 250.f; // Jack diff multiplier
    
    for (float jd : j) { // Iterate local jack difficulties
        // Decrease if jack difficulty is more than 133% player skill
        mod += ( (jd/x)*4/3 - 1) / 250.f;
        
        if (mod > 1.f)
            ceiling += (mod - 1) / fscale;
        
        // Clamp mod between 1 and 1.15
        mod = CalcClamp(mod, 1.f, base_ceiling * sqrt(ceiling));
        
        jd *= mod;
        
        if (x < jd) { // If player skill below jack diffiulty
            // This can cause output to decrease if 0.96 * i < x < i
            output += 1.f - pow(x / (jd * 0.96f), 1.5f);
        }
    }
    
    return CalcClamp(7.f * output, 0.f, 10000.f);
}

// Go through every note and determine a local jack speed difficulty at
// each place. That means:
//  1) taking the average of the last three note intervals,
//  2) (maybe bump that if the recent jack was really fast)
//  3) calculating 2800ms/interval_avg,
//  4) and maxing that out at the equivalent of 56 local NPS
// Returns a vector of each local jack speed difficulty
JackSeq Calc::SequenceJack(const vector<NoteInfo>& NoteInfo, unsigned int t, float music_rate) {
    vector<float> output;
    float last = -5.f;
    
    // Three most recent note intervals in ms. interval3 is the most
    // recent one.
    float interval1;
    float interval2 = 0.f;
    float interval3 = 0.f;
    
    unsigned int track = 1u << t;

    for (auto i : NoteInfo) {
        if (i.notes & track) { // If there's notes on the track
            float current_time = i.rowTime / music_rate;
            interval1 = interval2;
            interval2 = interval3;
            interval3 = 1000.f * (current_time - last);
            last = current_time;
            
            // Take the average of last three note intervals
            float interval_avg = (interval1 + interval2 + interval3) / 3.f;
            
            // If the last interval was really fast, jump to using that
            // instead of the average
            if (interval3 * 1.4 < interval_avg)
                interval_avg = interval3 * 1.4;
            
            // Difficulty for the 'local' jack speed
            // For example 1 NPS => 2.8; 2 NPS => 5.6; 10 NPS => 28
            float a = 2800.f / interval_avg;
            
            // Max out local jack speed difficulty at 56 NPS
            float v = min(a, 50.f);
            
            output.emplace_back(v);
        }
    }
    return output;
}

void Calc::InitializeHands(const vector<NoteInfo>& NoteInfo, float music_rate) {
    // Number of intervals
    numitv = static_cast<int>(std::ceil(NoteInfo.back().rowTime / (music_rate * IntervalSpan)));

    ProcessedFingers fingers;
    for (int i = 0; i < 4; i++) {
        fingers.emplace_back(ProcessFinger(NoteInfo, i, music_rate));
    }

    left_hand.InitDiff(fingers[0], fingers[1]);
    left_hand.InitPoints(fingers[0], fingers[1]);
    left_hand.ohjumpscale = OHJumpDownscaler(NoteInfo, 1, 2);
    left_hand.anchorscale = Anchorscaler(NoteInfo, 1, 2);
    left_hand.rollscale = RollDownscaler(fingers[0], fingers[1]);
    left_hand.hsscale = HSDownscaler(NoteInfo);
    left_hand.jumpscale = JumpDownscaler(NoteInfo);

    right_hand.InitDiff(fingers[2], fingers[3]);
    right_hand.InitPoints(fingers[2], fingers[3]);
    right_hand.ohjumpscale = OHJumpDownscaler(NoteInfo, 4, 8);
    right_hand.anchorscale = Anchorscaler(NoteInfo, 4, 8);
    right_hand.rollscale = RollDownscaler(fingers[2], fingers[3]);
    right_hand.hsscale = left_hand.hsscale;
    right_hand.jumpscale = left_hand.jumpscale;

    j0 = SequenceJack(NoteInfo, 0, music_rate);
    j1 = SequenceJack(NoteInfo, 1, music_rate);
    j2 = SequenceJack(NoteInfo, 2, music_rate);
    j3 = SequenceJack(NoteInfo, 3, music_rate);
}

Finger Calc::ProcessFinger(const vector<NoteInfo>& NoteInfo, unsigned int t, float music_rate) {
    int Interval = 0;
    float last = -5.f;
    Finger AllIntervals(numitv,vector<float>());
    if (t == 0)
        nervIntervals = vector<vector<int>>(numitv, vector<int>());
    unsigned int column = 1u << t;

    for (size_t i = 0; i < NoteInfo.size(); i++) {
        float scaledtime = NoteInfo[i].rowTime / music_rate;

        while (scaledtime > static_cast<float>(Interval + 1) * IntervalSpan)
            ++Interval;

        if (NoteInfo[i].notes & column) {
            AllIntervals[Interval].emplace_back(CalcClamp(1000 * (scaledtime - last), 40.f, 5000.f));
            last = scaledtime;
        }

        if (t == 0 && NoteInfo[i].notes != 0)
            nervIntervals[Interval].emplace_back(i);
    }
    return AllIntervals;
}

void Calc::TotalMaxPoints() {
    for (size_t i = 0; i < left_hand.v_itvpoints.size(); i++)
        MaxPoints += static_cast<float>(left_hand.v_itvpoints[i] + right_hand.v_itvpoints[i]);
}

float Calc::CalcScoreForPlayerSkill(float player_skill, ChiselFlags flags) {
    float gotpoints;
    if (flags & CHISEL_JACK) {
        // Max achievable points, minus the points the player's losing
        // from jack patterns
        gotpoints = MaxPoints
                - JackLoss(j0, player_skill)
                - JackLoss(j1, player_skill)
                - JackLoss(j2, player_skill)
                - JackLoss(j3, player_skill);
    } else {
        // Expected achieved points by left and right hand summed up
        gotpoints = left_hand.CalcInternal(player_skill, flags);
        gotpoints += right_hand.CalcInternal(player_skill, flags);
    }
    
    return gotpoints / MaxPoints;
}

// Approximate player skill required to achieve `score_goal`. The
// approximation can be influenced via the `flags`.
float Calc::Chisel(float player_skill, float resolution, float score_goal, ChiselFlags flags) {
    auto check_if_too_low = [this, flags, score_goal](float player_skill) {
        float score = CalcScoreForPlayerSkill(player_skill, flags);
        return score < score_goal;
    };
    return approximate(player_skill, resolution, 7, check_if_too_low, true);
}

// Looks at 6 smallest note intervals and returns 1375 / avg_interval_ms
// which could also be expressed as 1.375 * avg_nps.
float Hand::CalcMSEstimate(vector<float>& input) {
    if (input.empty())
        return 0.f;

    sort(input.begin(), input.end());
    float m = 0;
    input[0] *= 1.066f; //This is gross
    size_t End = min(input.size(), static_cast<size_t>(6));
    for (size_t i = 0; i < End; i++)
        m += input[i];
    return 1375.f * End / m;
}

void Hand::InitDiff(Finger& f1, Finger& f2) {
    v_itvNPSdiff = vector<float>(f1.size());
    v_itvMSdiff = vector<float>(f1.size());

    for (size_t i = 0; i < f1.size(); i++) {
        float nps = 1.6f * static_cast<float>(f1[i].size() + f2[i].size());
        float left_difficulty = CalcMSEstimate(f1[i]);
        float right_difficulty = CalcMSEstimate(f2[i]);
        float difficulty = max(left_difficulty, right_difficulty);
        v_itvNPSdiff[i] = finalscaler * nps;
        v_itvMSdiff[i] = finalscaler * (5.f * difficulty + 4.f * nps) / 9.f;
    }
    Smooth(v_itvNPSdiff, 0.f);
    if (SmoothDifficulty)
        DifficultyMSSmooth(v_itvMSdiff);
}

void Hand::InitPoints(const Finger& f1, const Finger& f2) {
    for (size_t i = 0; i < f1.size(); i++)
        v_itvpoints.emplace_back(static_cast<int>(f1[i].size()) + static_cast<int>(f2[i].size()));
}

void Hand::StamAdjust(float x, vector<float>& diff) {
    float floor = 1.f;          // stamina multiplier min (increases as chart advances)
    float mod = 1.f;            // multiplier
    float avs1;
    float avs2 = 0.f;

    for (float & i : diff) {
        avs1 = avs2;
        avs2 = i;
        float ebb = (avs1 + avs2) / 2;
        mod += ((ebb / (prop*x)) - 1) / mag;
        if (mod > 1.f)
            floor += (mod - 1) / fscale;
        mod = CalcClamp(mod, floor, ceil);
        i *= mod;
    }
}

// `nps`: Whether to use MS or NPS diff estimates
float Hand::CalcInternal(float x, ChiselFlags flags) {
    vector<float> diff = (flags & CHISEL_NPS) ? v_itvNPSdiff : v_itvMSdiff;
    
    for (size_t i = 0; i < diff.size(); ++i) {
        diff[i] *= anchorscale[i] * rollscale[i];
        
        if (flags & CHISEL_HS) {
            diff[i] *= sqrt(ohjumpscale[i]) * jumpscale[i];
        } else if (flags & CHISEL_JS) {
            diff[i] *= hsscale[i] * hsscale[i] * sqrt(ohjumpscale[i]) * jumpscale[i];
        } else if (flags & CHISEL_NPS) {
            diff[i] *= hsscale[i] * hsscale[i] * hsscale[i] * ohjumpscale[i] * ohjumpscale[i] * jumpscale[i] * jumpscale[i];
        } else {
            diff[i] *= sqrt(ohjumpscale[i]);
        }
    }

    if (flags & CHISEL_STAM)
        StamAdjust(x, diff);
    float output = 0.f;
    for (size_t i = 0; i < diff.size(); i++) {
        float delta = v_itvpoints[i];
        if (x <= diff[i])
            delta *= pow(x / diff[i], 1.8f);
        output += delta;
    }
    return output;
}

vector<float> Calc::OHJumpDownscaler(const vector<NoteInfo>& NoteInfo, unsigned int firstNote, unsigned int secondNote) {
    vector<float> output;

    std::cout << "sqrt(-4) = " << pow(-4, 0.5) << endl;

    for (const vector<int>& interval : nervIntervals) {
        int taps = 0;
        int jumps = 0;
        for (int row : interval) {
            int columns = 0;
            if (NoteInfo[row].notes & firstNote) {
                ++columns;
            }
            if (NoteInfo[row].notes & secondNote) {
                ++columns;
            }
            if (columns == 2) {
                jumps++;
                taps += 2; //this gets added twice intentionally to mimic mina's ratings more closely
            }
            taps += columns;
        }
        if (taps == 0) {
            output.push_back(1);
        } else {
            float jump_proportion = static_cast<float>(jumps) / static_cast<float>(taps);
            // When 62.5% of taps are jumps, the downscaler will reach 0
            output.push_back(pow(1 - (1.6f * jump_proportion), 0.25f));
        }

        if (logpatterns)
            std::cout << "ohj " << output.back() << std::endl;
    }

    if (SmoothPatterns)
        Smooth(output, 1.f);
    return output;
}

vector<float> Calc::Anchorscaler(const vector<NoteInfo>& NoteInfo, unsigned int firstNote, unsigned int secondNote) {
    vector<float> output(nervIntervals.size());

    for (size_t i = 0; i < nervIntervals.size(); i++) {
        int lcol = 0;
        int rcol = 0;
        for (int row : nervIntervals[i]) {
            if (NoteInfo[row].notes & firstNote)
                ++lcol;
            if (NoteInfo[row].notes & secondNote)
                ++rcol;
        }
        bool anyzero = lcol == 0 || rcol == 0;
        output[i] = anyzero ? 1.f : CalcClamp(sqrt(1 - (static_cast<float>(min(lcol, rcol)) / static_cast<float>(max(lcol, rcol)) / 4.45f)), 0.8f, 1.05f);

        fingerbias += (static_cast<float>(max(lcol, rcol)) + 2.f) / (static_cast<float>(min(lcol, rcol)) + 1.f);

        if (logpatterns)
            std::cout << "an " << output[i] << std::endl;
    }

    if (SmoothPatterns)
        Smooth(output, 1.f);
    return output;
}

vector<float> Calc::HSDownscaler(const vector<NoteInfo>& NoteInfo) {
    vector<float> output(nervIntervals.size());

    for (size_t i = 0; i < nervIntervals.size(); i++) {
        unsigned int taps = 0;
        unsigned int handtaps = 0;
        for (int row : nervIntervals[i]) {
            unsigned int notes = column_count(NoteInfo[row].notes);
            taps += notes;
            if (notes == 3)
                handtaps++;
        }
        output[i] = taps != 0 ? sqrt(sqrt(1 - (static_cast<float>(handtaps) / static_cast<float>(taps)))) : 1.f;

        if (logpatterns)
            std::cout << "hs " << output[i] << std::endl;
    }

    if (SmoothPatterns)
        Smooth(output, 1.f);
    return output;
}

vector<float> Calc::JumpDownscaler(const vector<NoteInfo>& NoteInfo) {
    vector<float> output(nervIntervals.size());

    for (size_t i = 0; i < nervIntervals.size(); i++) {
        unsigned int taps = 0;
        unsigned int jumps = 0;
        for (int row : nervIntervals[i]) {
            unsigned int notes = column_count(NoteInfo[row].notes);
            taps += notes;
            if (notes == 2)
                jumps++;
        }
        output[i] = taps != 0 ? sqrt(sqrt(1 - (static_cast<float>(jumps) / static_cast<float>(taps) / 3.f))) : 1.f;

        if (logpatterns)
            std::cout << "ju " << output[i] << std::endl;
    }
    if (SmoothPatterns)
        Smooth(output, 1.f);
    return output;
}


vector<float> Calc::RollDownscaler(const Finger& f1, const Finger& f2) {
    vector<float> output(f1.size());    //this is slightly problematic because if one finger is longer than the other
                                        //you could potentially have different results with f1 and f2 switched
    for (size_t i = 0; i < f1.size(); i++) {
        if (f1[i].size() + f2[i].size() <= 1) {
            output[i] = 1.f;
            continue;
        }
        vector<float> hand_intervals;
        for (float time1 : f1[i])
            hand_intervals.emplace_back(time1);
        for (float time2 : f2[i])
            hand_intervals.emplace_back(time2);

        float interval_mean = mean(hand_intervals);

        for (float & note : hand_intervals)
            if (interval_mean / note < 0.6f)
                note = interval_mean;

        float interval_cv = cv(hand_intervals) + 0.85f;
        output[i] = interval_cv >= 1.0f ? min(sqrt(sqrt(interval_cv)), 1.075f) : interval_cv*interval_cv*interval_cv;

        if (logpatterns)
            std::cout << "ro " << output[i] << std::endl;
    }

    if (SmoothPatterns)
        Smooth(output, 1.f);

    return output;
}

// Function to generate SSR rating
DifficultyRating MinaSDCalc(const vector<NoteInfo>& NoteInfo, float musicrate, float goal) {
    if (NoteInfo.empty()) {
        return DifficultyRating {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    }
    return std::make_unique<Calc>()->CalcMain(NoteInfo, musicrate, goal);
}

// Wrap difficulty calculation for all rates from 0.7 to 2.1, with 0.1
// step
MinaSD MinaSDCalc(const vector<NoteInfo>& NoteInfo) {
    MinaSD allrates;
    int lower_rate = 7;
    int upper_rate = 21;

    if (!NoteInfo.empty())
        for (int i = lower_rate; i < upper_rate; i++)
            allrates.emplace_back(MinaSDCalc(NoteInfo, static_cast<float>(i) / 10.f, 0.93f));
    else {
        DifficultyRating output{0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
        for (int i = lower_rate; i < upper_rate; i++)
            allrates.emplace_back(output);
    }
    return allrates;
}

int GetCalcVersion()
{
    return -1;
}
