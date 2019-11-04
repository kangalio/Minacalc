// MinaCalc.cpp : Defines the exported functions for the DLL application.
//

#include "minacalc.h"
#include <iostream>
#include <algorithm>
#include <thread>
#include <mutex>
#include <cmath>


#define SAFE_DELETE(p){ delete p; p = NULL;}

template<typename T>
T CalcClamp(T x, T l, T h) {
    if (x > h)
        x = h;
    else
    if (x < l)
        x = l;
return x;
}

inline float mean(vector<float>& input) {
    float sum = 0.f;
    for (float i : input)
        sum += i;

    return sum / input.size();
}

// Coefficient of variation
inline float cv(vector<float> &input) {
    float sd = 0.f;
    float average = mean(input);
    for (float i : input)
        sd += (i - average)*(i - average);

    return sqrt(sd / input.size()) / average;
}

inline float downscale_low_accuracy_scores(float f, float sg) {
    if (sg >= 0.93f)
        return f;
    return min(max(f - sqrt(0.93f - sg), 0.f), 100.f);
}

// Specifically for pattern modifiers as the neutral value is 1
inline void PatternSmooth(vector<float>& input) {
    float f1 = 1.f;
    float f2 = 1.f;
    float f3 = 1.f;
    float total = 3.f;

    for (float & i : input) {
        total -= f1;
        f1 = f2;
        f2 = f3;
        f3 = i;
        total += f3;
        i = (f1 + f2 + f3) / 3;
    }
}

inline void DifficultySmooth(vector<float>& input) {
    float f1 = 0.f;
    float f2 = 0.f;
    float f3 = 0.f;
    float total = 0.f;

    for (float & i : input) {
        total -= f1;
        f1 = f2;
        f2 = f3;
        f3 = i;
        total += f3;
        i = (f1 + f2 + f3) / 3;
    }
}

inline void DifficultyMSSmooth(vector<float>& input) {
    float f1;
    float f2 = 0.f;

    for (float & i : input) {
        f1 = f2;
        f2 = i;
        i = (f1 + f2) / 2.f;
    }
}

inline float AggregateScores(vector<float>& invector, float rating, float res, int iter) {
    float sum;
    do {
        rating += res;
        sum = 0.0f;
        for (float i : invector) {
            sum += 2.f / erfc(0.5f*(i - rating)) - 1.f;
        }
    } while (3 < sum);
    if (iter == 11)
        return rating;
    return AggregateScores(invector, rating - res, res / 2.f, iter + 1);
}

float normalizer(float x, float y, float z1, float z2) {
    float norm = CalcClamp(((x / y) - 1.f)*z1, 0.f, 1.f);
    return x * z2 * norm + x * (1.f - z2);
}

float jumpprop(const vector<NoteInfo>& NoteInfo) {
    int left = 1;
    int down = 2;
    int up = 4;
    int right = 8;

    int taps = 0;
    int jumps = 0;

    for (auto r : NoteInfo) {
        int notes = (r.notes & left ? 1 : 0) + (r.notes & down ? 1 : 0) + (r.notes & up ? 1 : 0) + (r.notes & right ? 1 : 0);
        taps += notes;
        if (notes == 2)
            jumps += notes;
    }

    return static_cast<float>(jumps) / static_cast<float>(taps);
}

float handprop(const vector<NoteInfo>& NoteInfo) {
    int left = 1;
    int down = 2;
    int up = 4;
    int right = 8;

    int taps = 0;
    int hands = 0;

    for (auto r : NoteInfo) {
        int notes = (r.notes & left ? 1 : 0) + (r.notes & down ? 1 : 0) + (r.notes & up ? 1 : 0) + (r.notes & right ? 1 : 0);
        taps += notes;
        if (notes == 3)
            hands += notes;
    }

    return static_cast<float>(hands) / static_cast<float>(taps);
}

float quadprop(const vector<NoteInfo>& NoteInfo) {
    int left = 1;
    int down = 2;
    int up = 4;
    int right = 8;

    int taps = 0;
    int quads = 0;

    for (auto r : NoteInfo) {
        int notes = (r.notes & left ? 1 : 0) + (r.notes & down ? 1 : 0) + (r.notes & up ? 1 : 0) + (r.notes & right ? 1 : 0);
        taps += notes;
        if (notes == 4)
            quads += notes;
    }

    return static_cast<float>(quads) / static_cast<float>(taps);
}

vector<float> skillset_vector(DifficultyRating& difficulty) {
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
    return max(difficulty.stream,
            max(difficulty.jumpstream,
                    max(difficulty.handstream,
                            max( difficulty.stamina,
                                    max(difficulty.jack,
                                            max(difficulty.chordjack, difficulty.technical))))));
}

DifficultyRating Calc::CalcMain(const vector<NoteInfo>& NoteInfo) {
    float grindscaler = CalcClamp(0.93f + (0.07f * (NoteInfo.back().rowTime - 30.f) / 30.f), 0.93f, 1.f)
            * CalcClamp(0.873f + (0.13f * (NoteInfo.back().rowTime - 15.f) / 15.f), 0.87f, 1.f);

    float shortstamdownscaler = CalcClamp(0.9f + (0.1f * (NoteInfo.back().rowTime - 150.f) / 150.f), 0.9f, 1.f);

    float jprop = jumpprop(NoteInfo);
    float nojumpsdownscaler = CalcClamp(0.8f + (0.2f * (jprop + 0.5f)), 0.8f, 1.f);
    float manyjumpsdownscaler = CalcClamp(1.43f - jprop, 0.85f, 1.f);

    float hprop = handprop(NoteInfo);
    float nohandsdownscaler = CalcClamp(0.8f + (0.2f * (hprop + 0.75f)), 0.8f, 1.f);
    float allhandsdownscaler = CalcClamp(1.23f - hprop, 0.85f, 1.f);

    float qprop = quadprop(NoteInfo);
    float lotquaddownscaler = CalcClamp(1.13f - qprop, 0.85f, 1.f);

    float jumpthrill = CalcClamp(1.625f - jprop - hprop, 0.85f, 1.f);

    InitializeHands(NoteInfo);
    TotalMaxPoints();
    float stream = Chisel(0.1f, 10.24f, false, false, true, false, false);
    float js = Chisel(0.1f, 10.24f, false, false, true, true, false);
    float hs = Chisel(0.1f, 10.24f, false, false, true, false, true);
    float tech = Chisel(0.1f, 10.24f, false, false, false, false, false);
    float jack = Chisel(0.1f, 10.24f, false, true, true, false, false);

    float techbase = max(stream, jack);
    tech = CalcClamp((tech / techbase)*tech, tech * 0.85f, tech);

    float stam;
    if (stream > tech || js > tech || hs > tech)
        if (stream > js && stream > hs)
            stam = Chisel(stream - 0.1f, 2.56f, true, false, true, false, false);
        else if (js > hs)
            stam = Chisel(js - 0.1f, 2.56f, true, false, true, true, false);
        else
            stam = Chisel(hs - 0.1f, 2.56f, true, false, true, false, true);
    else
        stam = Chisel(tech - 0.1f, 2.56f, true, false, false, false, false);

    js = normalizer(js, stream, 7.25f, 0.25f);
    hs = normalizer(hs, stream, 6.5f, 0.3f);
    hs = normalizer(hs, js, 11.5f, 0.15f);

    float stambase = max(max(stream, tech* 0.96f), max(js, hs));
    if (stambase == stream)
        stambase *= 0.975f;

    stam = normalizer(stam, stambase, 7.75f, 0.2f);

    float chordjack = jack * 0.75f;
    float technorm = max(max(stream, js), hs);
    tech = normalizer(tech, technorm, 8.f, .15f) * techscaler;

    DifficultyRating difficulty = DifficultyRating {0.0,
                                                    downscale_low_accuracy_scores(stream, Scoregoal),
                                                    downscale_low_accuracy_scores(js, Scoregoal),
                                                    downscale_low_accuracy_scores(hs, Scoregoal),
                                                    downscale_low_accuracy_scores(stam, Scoregoal),
                                                    downscale_low_accuracy_scores(jack, Scoregoal),
                                                    downscale_low_accuracy_scores(chordjack, Scoregoal),
                                                    downscale_low_accuracy_scores(tech, Scoregoal)
    };

    // chordjack
    float cj = difficulty.handstream;

    difficulty.stream *= allhandsdownscaler * manyjumpsdownscaler * lotquaddownscaler;
    difficulty.jumpstream *= nojumpsdownscaler * allhandsdownscaler * lotquaddownscaler;
    difficulty.handstream *= nohandsdownscaler * allhandsdownscaler * 1.015f * manyjumpsdownscaler * lotquaddownscaler;
    difficulty.stamina = CalcClamp(difficulty.stamina * shortstamdownscaler * 0.985f * lotquaddownscaler, 1.f,
                                   max(max(difficulty.stream, difficulty.jack), max(difficulty.jumpstream, difficulty.handstream)) * 1.1f);
    difficulty.technical *= allhandsdownscaler * manyjumpsdownscaler * lotquaddownscaler * 1.01f;

    cj = normalizer(cj, difficulty.handstream, 5.5f, 0.3f) * CalcClamp(qprop + hprop + jprop + 0.2f, 0.5f, 1.f) * 1.025f;

    if (cj > difficulty.jack)
        difficulty.chordjack = cj;
    else
        difficulty.chordjack *= 0.9f;


    dumbvalue = (dumbvalue / static_cast<float>(dumbcounter));
    float stupidvalue = CalcClamp(1.f - (dumbvalue - 2.55f), 0.85f, 1.f);
    difficulty.technical *= stupidvalue;

    if (stupidvalue <= 0.95f) {
        difficulty.jack *= 1.f + (1.f - sqrt(stupidvalue));
    }

    float skadoot = max(difficulty.handstream, difficulty.jumpstream);
    if (difficulty.stream < skadoot)
        difficulty.stream -= sqrt(skadoot - difficulty.stream);

    vector<float> temp_vec = skillset_vector(difficulty);
    float overall = AggregateScores(temp_vec, 0.f, 10.24f, 1);
    difficulty.overall = downscale_low_accuracy_scores(overall, Scoregoal);

    temp_vec = skillset_vector(difficulty);
    float aDvg = mean(temp_vec) * 1.2f;
    difficulty.overall = downscale_low_accuracy_scores(min(difficulty.overall, aDvg) * grindscaler, Scoregoal);
    difficulty.stream = downscale_low_accuracy_scores(min(difficulty.stream, aDvg * 1.0416f) * grindscaler, Scoregoal);
    difficulty.jumpstream = downscale_low_accuracy_scores(min(difficulty.jumpstream, aDvg * 1.0416f) * grindscaler, Scoregoal) * jumpthrill;
    difficulty.handstream = downscale_low_accuracy_scores(min(difficulty.handstream, aDvg) * grindscaler, Scoregoal) * jumpthrill;
    difficulty.stamina = downscale_low_accuracy_scores(min(difficulty.stamina, aDvg) * grindscaler, Scoregoal) * sqrt(jumpthrill) * 0.996f;
    difficulty.jack = downscale_low_accuracy_scores(min(difficulty.jack, aDvg) * grindscaler, Scoregoal);
    difficulty.chordjack = downscale_low_accuracy_scores(min(difficulty.chordjack, aDvg) * grindscaler, Scoregoal);
    difficulty.technical = downscale_low_accuracy_scores(min(difficulty.technical, aDvg * 1.0416f) * grindscaler, Scoregoal) * sqrt(jumpthrill);

    float highest = max(difficulty.overall, highest_difficulty(difficulty));

    vector<float> temp = skillset_vector(difficulty);
    difficulty.overall = AggregateScores(temp, 0.f, 10.24f, 1);

    float dating = CalcClamp(0.5f + (highest / 100.f), 0.f, 0.9f);

    if (Scoregoal < dating) {
        difficulty = DifficultyRating {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
    }

    difficulty.jack *= 1.0075f;

    if (highest == difficulty.technical) {
        difficulty.technical -= CalcClamp(4.5f - difficulty.technical + difficulty.handstream, 0.f, 4.5f);
        difficulty.technical -= CalcClamp(4.5f - difficulty.technical + difficulty.jumpstream, 0.f, 4.5f);
    }

    difficulty.technical *= 1.025f;
    difficulty.overall = highest_difficulty(difficulty);

    return difficulty;
}

// ugly jack stuff
vector<float> Calc::JackStamAdjust(vector<float>& j, float x) {
    vector<float> output(j.size());
    float floor = 1.f;
    float mod = 1.f;
    float ceil = 1.15f;
    float fscale = 1750.f;
    float prop = 0.75f;
    float mag = 250.f;
    float multstam = 1.f;

    for (size_t i = 0; i < j.size(); i++) {
        mod += ((j[i] * multstam / (prop*x)) - 1) / mag;
        if (mod > 1.f)
            floor += (mod - 1) / fscale;
        mod = CalcClamp(mod, 1.f, ceil * sqrt(floor));
        output[i] = j[i] * mod;
    }
    return output;
}

float Calc::JackLoss(vector<float>& j, float x) {
    const vector<float>& v = JackStamAdjust(j, x);
    float output = 0.f;
    for (float i : v) {
        if (x < i)
            output += 7.f - (7.f * pow(x / (i * 0.96f), 1.5f));
    }
    return CalcClamp(output, 0.f, 10000.f);
}

JackSeq Calc::SequenceJack(const vector<NoteInfo>& NoteInfo, int t) {
    vector<float> output;
    float last = -5.f;
    float mats1;
    float mats2 = 0.f;
    float mats3 = 0.f;
    float timestamp;
    int track = 1 << t;

    for (auto i : NoteInfo) {
        float scaledtime = i.rowTime / MusicRate;
        if (i.notes & track) {
            mats1 = mats2;
            mats2 = mats3;
            mats3 = 1000.f * (scaledtime - last);
            last = scaledtime;
            timestamp = CalcClamp((mats1 + mats2 + mats3) / 3.f, 25.f, mats3 * 1.4f);
            output.emplace_back(CalcClamp(1 / timestamp * 2800.f, 0.f, 50.f));
        }
    }
    return output;
}

int Calc::fastwalk(const vector<NoteInfo>& NoteInfo) {
    int Interval = 0;
    for (auto i : NoteInfo) {
        if (i.rowTime / MusicRate >= Interval * IntervalSpan)
            ++Interval;
    }
    return Interval;
}


void Calc::InitializeHands(const vector<NoteInfo>& NoteInfo) {
    numitv = fastwalk(NoteInfo);

    ProcessedFingers left_fingers;
    ProcessedFingers right_fingers;
    for (int i = 0; i < 4; i++) {
        if (i < 2)
            left_fingers.emplace_back(ProcessFinger(NoteInfo, i));
        else
            right_fingers.emplace_back(ProcessFinger(NoteInfo, i));
    }

    left_hand = new Hand;
    left_hand->InitHand(left_fingers[0], left_fingers[1]);
    left_hand->ohjumpscale = OHJumpDownscaler(NoteInfo, 1, 2);
    left_hand->anchorscale = Anchorscaler(NoteInfo, 1, 2);
    left_hand->rollscale = RollDownscaler(left_fingers[0], left_fingers[1]);
    left_hand->hsscale = HSDownscaler(NoteInfo);
    left_hand->jumpscale = JumpDownscaler(NoteInfo);

    right_hand = new Hand;
    right_hand->InitHand(right_fingers[0], right_fingers[1]);
    right_hand->ohjumpscale = OHJumpDownscaler(NoteInfo, 4, 8);
    right_hand->anchorscale = Anchorscaler(NoteInfo, 4, 8);
    right_hand->rollscale = RollDownscaler(right_fingers[0], right_fingers[1]);
    right_hand->hsscale = left_hand->hsscale;
    right_hand->jumpscale = left_hand->jumpscale;

    j0 = SequenceJack(NoteInfo, 0);
    j1 = SequenceJack(NoteInfo, 1);
    j2 = SequenceJack(NoteInfo, 2);
    j3 = SequenceJack(NoteInfo, 3);

    vector<Finger> ltmp;
    vector<Finger> rtmp;
    left_fingers.swap(ltmp);
    right_fingers.swap(rtmp);

    left_fingers.shrink_to_fit();
    right_fingers.shrink_to_fit();
}

Finger Calc::ProcessFinger(const vector<NoteInfo>& NoteInfo, int t) {
    int Interval = 1;
    float last = -5.f;
    Finger AllIntervals(numitv);
    vector<float> CurrentInterval;
    vector<int> itvnervtmp;
    vector<vector<int>> itvnerv(numitv);

    int left = 1;
    int down = 2;
    int up = 4;
    int right = 8;

    int column = 1 << t;
    for (size_t i = 0; i < NoteInfo.size(); i++) {
        float scaledtime = NoteInfo[i].rowTime / MusicRate;

        if (scaledtime >= Interval * IntervalSpan) {
            AllIntervals[Interval - 1] = CurrentInterval;
            CurrentInterval.clear();

            itvnerv[Interval - 1] = itvnervtmp;
            itvnervtmp.clear();
            ++Interval;
        }

        if (NoteInfo[i].notes & column) {
            CurrentInterval.emplace_back(CalcClamp(1000 * (scaledtime - last), 40.f, 5000.f));
            last = scaledtime;
        }

        if (t == 0 && (NoteInfo[i].notes & left || NoteInfo[i].notes & down || NoteInfo[i].notes & up || NoteInfo[i].notes & right))
        {
            itvnervtmp.emplace_back(i);
        }
    }

    if(t == 0)
        nervIntervals = itvnerv;
    return AllIntervals;
}

void Calc::TotalMaxPoints() {
    for (size_t i = 0; i < left_hand->v_itvpoints.size(); i++)
        MaxPoints += static_cast<int>(left_hand->v_itvpoints[i] + right_hand->v_itvpoints[i]);
}

float Calc::Chisel(float player_skill, float resolution, bool stamina, bool jack, bool nps, bool js, bool hs) {
    float gotpoints;
    for (int iter = 1; iter <= 7; iter++) {
        do {
            if (player_skill > 100.f)
                return player_skill;
            player_skill += resolution;
            if (jack) {
                gotpoints = MaxPoints - JackLoss(j0, player_skill) - JackLoss(j1, player_skill) - JackLoss(j2, player_skill) -
                            JackLoss(j3, player_skill);
            } else
                gotpoints = left_hand->CalcInternal(player_skill, stamina, nps, js, hs) +
                            right_hand->CalcInternal(player_skill, stamina, nps, js, hs);

        } while (gotpoints / MaxPoints < Scoregoal);
        player_skill -= resolution;
        resolution /= 2.f;
    }
    return player_skill + 2.f * resolution;
}


// Hand stuff
void Hand::InitHand(Finger & f1, Finger & f2) {
    InitDiff(f1, f2);
    InitPoints(f1, f2);
}

float Hand::CalcMSEstimate(vector<float>& input) {
    if (input.empty())
        return 0.f;

    sort(input.begin(), input.end());
    float m = 0;
    input[0] *= 1.066f;
    size_t End = min(input.size(), static_cast<size_t>(6));
    for (size_t i = 0; i < End; i++)
        m += input[i];
    return 1 / (m / (End)) * 1375;
}

void Hand::InitDiff(Finger& f1, Finger& f2) {
    vector<float> tmpNPS(f1.size());
    vector<float> tmpMS(f1.size());

    for (size_t i = 0; i < f1.size(); i++) {
        float nps = 1.6f * static_cast<float>(f1[i].size() + f2[i].size());		// intervalspan
        float aa = CalcMSEstimate(f1[i]);
        float bb = CalcMSEstimate(f2[i]);
        float ms = max(aa, bb);
        tmpNPS[i] = finalscaler * nps;
        tmpMS[i] = finalscaler * (5.f * ms + 4.f * nps) / 9.f;
    }
    if (SmoothDifficulty)
        DifficultyMSSmooth(tmpMS);

    DifficultySmooth(tmpNPS);
    v_itvNPSdiff = tmpNPS;
    v_itvMSdiff = tmpMS;
}

void Hand::InitPoints(Finger& f1, Finger& f2) {
    for (size_t i = 0; i < f1.size(); i++)
        v_itvpoints.emplace_back(static_cast<int>(f1[i].size()) + static_cast<int>(f2[i].size()));
}

vector<float> Hand::StamAdjust(float x, vector<float> diff) {
    vector<float> output(diff.size());
    float floor = 1.f;			// stamina multiplier min (increases as chart advances)
    float mod = 1.f;			// mutliplier

    float avs1;
    float avs2 = 0.f;

    for (size_t i = 0; i < diff.size(); i++) {
        avs1 = avs2;
        avs2 = diff[i];
        float ebb = (avs1 + avs2) / 2;
        mod += ((ebb / (prop*x)) - 1) / mag;
        if (mod > 1.f)
            floor += (mod - 1) / fscale;
        mod = CalcClamp(mod, floor, ceil);
        output[i] = diff[i] * mod;
    }
    return output;
}

float Hand::CalcInternal(float x, bool stam, bool nps, bool js, bool hs) {
    vector<float> diff;

    if (nps)
        diff = v_itvNPSdiff;
    else
        diff = v_itvMSdiff;

    for (size_t i = 0; i < diff.size(); ++i) {
        if (hs)
            diff[i] = diff[i] * anchorscale[i] * sqrt(ohjumpscale[i]) * rollscale[i] * jumpscale[i];
        else if (js)
            diff[i] = diff[i] * pow(hsscale[i], 2) * anchorscale[i] * sqrt(ohjumpscale[i]) * rollscale[i] * jumpscale[i];
        else if (nps)
            diff[i] = diff[i] * pow(hsscale[i], 3) * anchorscale[i] * pow(ohjumpscale[i], 2) * rollscale[i] * pow(jumpscale[i], 2);
        else
            diff[i] = diff[i] * anchorscale[i] * sqrt(ohjumpscale[i]) * rollscale[i];
    }

    const vector<float>& v = stam ? StamAdjust(x, diff) : diff;
    //dum = v;
    float output = 0.f;
    for (size_t i = 0; i < v.size(); i++) {
        if (x > v[i])
            output += v_itvpoints[i];
        else
            output += v_itvpoints[i] * pow(x / v[i], 1.8f);
    }
    return output;
}


// pattern modifiers
vector<float> Calc::OHJumpDownscaler(const vector<NoteInfo>& NoteInfo, int firstNote, int secondNote) {
    vector<float> output(nervIntervals.size());

    for (size_t i = 0; i < nervIntervals.size(); i++) {
        int taps = 0;
        int jumptaps = 0;
        for (int row : nervIntervals[i]) {
            if (NoteInfo[row].notes & firstNote) {
                ++taps;
                if (NoteInfo[row].notes & secondNote) {
                    jumptaps += 2;
                    ++taps;
                }
            }
        }
        output[i] = taps != 0 ? pow(1 - (static_cast<float>(jumptaps) / static_cast<float>(taps) / 2.5f), 0.25f) : 1.f;

        if (logpatterns)
            cout << "ohj " << output[i] << endl;
    }

    if (SmoothPatterns)
        PatternSmooth(output);
    return output;
}

// pattern modifiers
vector<float> Calc::Anchorscaler(const vector<NoteInfo>& NoteInfo, int firstNote, int secondNote) {
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

        float stupidthing = (static_cast<float>(max(lcol, rcol)) + 2.f) / (static_cast<float>(min(lcol, rcol)) + 1.f);
        dumbvalue += stupidthing;
        ++dumbcounter;

        if (logpatterns)
            cout << "an " << output[i] << endl;
    }

    if (SmoothPatterns)
        PatternSmooth(output);
    return output;
}


vector<float> Calc::HSDownscaler(const vector<NoteInfo>& NoteInfo) {
    vector<float> output(nervIntervals.size());
    int left = 1;
    int down = 2;
    int up = 4;
    int right = 8;

    for (size_t i = 0; i < nervIntervals.size(); i++) {
        if (nervIntervals[i].empty())
            output[i] = 1.f;
        else {
            int taps = 0;
            int handtaps = 0;
            for (int row : nervIntervals[i]) {
                int notes = (NoteInfo[row].notes & left ? 1 : 0) + (NoteInfo[row].notes & down ? 1 : 0) + (NoteInfo[row].notes & up ? 1 : 0) + (NoteInfo[row].notes & right ? 1 : 0);
                taps += notes;
                if (notes == 3)
                    handtaps += notes;
            }
            output[i] = taps != 0 ? sqrt(sqrt(1 - (static_cast<float>(handtaps) / static_cast<float>(taps) / 3.f))) : 1.f;

            if (logpatterns)
                cout << "hs " << output[i] << endl;
        }
    }

    if (SmoothPatterns)
        PatternSmooth(output);
    return output;
}

vector<float> Calc::JumpDownscaler(const vector<NoteInfo>& NoteInfo) {
    vector<float> output(nervIntervals.size());
    int left = 1;
    int down = 2;
    int up = 4;
    int right = 8;

    for (size_t i = 0; i < nervIntervals.size(); i++) {
        if (nervIntervals[i].empty())
            output[i] = 1.f;
        else {
            int taps = 0;
            int jumps = 0;
            for (int row : nervIntervals[i]) {
                int notes = (NoteInfo[row].notes & left ? 1 : 0) + (NoteInfo[row].notes & down ? 1 : 0) + (NoteInfo[row].notes & up ? 1 : 0) + (NoteInfo[row].notes & right ? 1 : 0);
                taps += notes;
                if (notes == 2)
                    jumps += notes;
            }
            output[i] = taps != 0 ? sqrt(sqrt(1 - (static_cast<float>(jumps) / static_cast<float>(taps) / 6.f))) : 1.f;

            if (logpatterns)
                cout << "ju " << output[i] << endl;
        }
    }
    if (SmoothPatterns)
        PatternSmooth(output);
    return output;
}


vector<float> Calc::RollDownscaler(Finger f1, Finger f2) {
    vector<float> output(f1.size());
    for (size_t i = 0; i < f1.size(); i++) {
        if (f1[i].empty() && f2[i].empty())
            output[i] = 1.f;
        else {
            vector<float> cvint;
            for (float & ii : f1[i])
                cvint.emplace_back(ii);
            for (float & ii : f2[i])
                cvint.emplace_back(ii);

            float mmm = mean(cvint);

            for (float & i : cvint)
                i = mmm / i < 0.6f ? mmm : i;


            if (cvint.size() == 1) {
                output[i] = 1.f;
                continue;
            }

            float dacv = cv(cvint);
            if (dacv >= 0.15)
                output[i] = sqrt(sqrt(0.85f + dacv));
            else
                output[i] = pow(0.85f + dacv, 3);
            output[i] = CalcClamp(output[i], 0.f, 1.075f);

            if (logpatterns)
                cout << "ro " << output[i] << endl;
        }
    }

    if (SmoothPatterns)
        PatternSmooth(output);

    return output;
}


void Calc::Purge() {
    vector<float> tmp1;
    vector<float> tmp2;
    vector<float> tmp3;
    vector<float> tmp4;

    j0.swap(tmp1);
    j1.swap(tmp2);
    j2.swap(tmp3);
    j3.swap(tmp4);

    j0.shrink_to_fit();
    j1.shrink_to_fit();
    j2.shrink_to_fit();
    j3.shrink_to_fit();

    vector<float> l1;
    vector<float> l2;
    vector<float> l3;
    vector<float> l4;
    vector<float> l5;

    left_hand->ohjumpscale.swap(l1);
    left_hand->anchorscale.swap(l2);
    left_hand->rollscale.swap(l3);
    left_hand->hsscale.swap(l4);
    left_hand->jumpscale.swap(l5);

    left_hand->ohjumpscale.shrink_to_fit();
    left_hand->anchorscale.shrink_to_fit();
    left_hand->rollscale.shrink_to_fit();
    left_hand->hsscale.shrink_to_fit();
    left_hand->jumpscale.shrink_to_fit();

    vector<float> r1;
    vector<float> r2;
    vector<float> r3;
    vector<float> r4;
    vector<float> r5;

    right_hand->ohjumpscale.swap(l1);
    right_hand->anchorscale.swap(l2);
    right_hand->rollscale.swap(l3);
    right_hand->hsscale.swap(l4);
    right_hand->jumpscale.swap(l5);

    right_hand->ohjumpscale.shrink_to_fit();
    right_hand->anchorscale.shrink_to_fit();
    right_hand->rollscale.shrink_to_fit();
    right_hand->hsscale.shrink_to_fit();
    right_hand->jumpscale.shrink_to_fit();


    SAFE_DELETE(left_hand);
    SAFE_DELETE(right_hand);
}

// Function to generate SSR rating
DifficultyRating MinaSDCalc(const vector<NoteInfo>& NoteInfo, float musicrate, float goal) {
    unique_ptr<Calc> doot = make_unique<Calc>();
    doot->MusicRate = musicrate;
    goal = CalcClamp(goal, 0.f, 0.965f);	// cap SSR at 96% so things don't get out of hand
    doot->Scoregoal = goal;
    DifficultyRating output = doot->CalcMain(NoteInfo);

    doot->Purge();

    return output;
}

// Wrap difficulty calculation for all standard rates
MinaSD MinaSDCalc(const vector<NoteInfo>& NoteInfo) {
    MinaSD allrates;
    int lower_rate = 7;
    int upper_rate = 21;

    if (!NoteInfo.empty())
        for (int i = lower_rate; i < upper_rate; i++)
            allrates.emplace_back(MinaSDCalc(NoteInfo,i / 10.f, 0.93f));
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