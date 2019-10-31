// MinaCalc.cpp : Defines the exported functions for the DLL application.
//

#include "minacalc.h"
#include <iostream>
#include <algorithm>
#include <thread>
#include <mutex>
#include <cmath>

using namespace std;


#define SAFE_DELETE(p){ delete p; p = NULL;}

template<typename T, typename U, typename V>
inline void CalcClamp(T& x, U l, V h) {
    if (x > static_cast<T>(h))
        x = static_cast<T>(h);
    else
    if (x < static_cast<T>(l))
        x = static_cast<T>(l);
}

inline float mean(vector<float>& input) {
    float sum = 0.f;
    for (float i : input)
        sum += i;

    return sum / input.size();
}

// Coefficient of variance
inline float cv(vector<float> &input) {
    float sum = 0.f;
    float mean;
    float sd = 0.f;

    for (float i : input)
        sum += i;

    mean = sum / input.size();
    for (float i : input)
        sd += pow(i - mean, 2);

    return sqrt(sd / input.size()) / mean;
}

inline float downscale_low_accuracy_scores(float& f, float sg) {
    CalcClamp(f, 0.f, 100.f);
    if (sg >= 0.93f)
        return f;
    float output = f * 1 - sqrt(0.93f - sg);
    CalcClamp(f, 0.f, 100.f);
    return output;
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
    float norm = ((x / y) - 1.f)*z1;
    CalcClamp(norm, 0.f, 1.f);
    float output = x * z2 * norm + x * (1.f - z2);
    return output;
}

float jumpprop(const vector<NoteInfo>& NoteInfo) {
    int left = 1;
    int down = 1 << 1;
    int up = 1 << 2;
    int right = 1 << 3;

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
    int down = 1 << 1;
    int up = 1 << 2;
    int right = 1 << 3;

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
    int down = 1 << 1;
    int up = 1 << 2;
    int right = 1 << 3;

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

float highest_difficulty(const DifficultyRating& difficulty) {
    return max(difficulty.stream,
            max(difficulty.jumpstream,
                    max(difficulty.handstream,
                            max( difficulty.stamina,
                                    max(difficulty.jack,
                                            max(difficulty.chordjack, difficulty.technical))))));
}

DifficultyRating Calc::CalcMain(const vector<NoteInfo>& NoteInfo) {
    //LOG->Trace("%f", etaner.back());
    float grindscaler = 0.93f + (0.07f * (NoteInfo.back().rowTime - 30.f) / 30.f);
    CalcClamp(grindscaler, 0.93f, 1.f);

    float grindscaler2 = 0.873f + (0.13f * (NoteInfo.back().rowTime - 15.f) / 15.f);
    CalcClamp(grindscaler2, 0.87f, 1.f);

    float shortstamdownscaler = 0.9f + (0.1f * (NoteInfo.back().rowTime - 150.f) / 150.f);
    CalcClamp(shortstamdownscaler, 0.9f, 1.f);

    float jprop = jumpprop(NoteInfo);
    float nojumpsdownscaler = 0.8f + (0.2f * (jprop + 0.5f));
    CalcClamp(nojumpsdownscaler, 0.8f, 1.f);

    float hprop = handprop(NoteInfo);

    float nohandsdownscaler = 0.8f + (0.2f * (hprop + 0.75f));
    CalcClamp(nohandsdownscaler, 0.8f, 1.f);

    float allhandsdownscaler = 1.23f - hprop;
    CalcClamp(allhandsdownscaler, 0.85f, 1.f);

    float manyjumpsdownscaler = 1.43f - jprop;
    CalcClamp(manyjumpsdownscaler, 0.85f, 1.f);

    float qprop = quadprop(NoteInfo);
    float lotquaddownscaler = 1.13f - qprop;
    CalcClamp(lotquaddownscaler, 0.85f, 1.f);

    float jumpthrill = 1.625f - jprop - hprop;
    CalcClamp(jumpthrill, 0.85f, 1.f);

    vector<float> output;
    output.reserve(8);

    InitializeHands(NoteInfo);
    TotalMaxPoints();
    float stream = Chisel(0.1f, 10.24f, 1, false, false, true, false, false);
    float js = Chisel(0.1f, 10.24f, 1, false, false, true, true, false);
    float hs = Chisel(0.1f, 10.24f, 1, false, false, true, false, true);
    float tech = Chisel(0.1f, 10.24f, 1, false, false, false, false, false);
    float jack = Chisel(0.1f, 10.24f, 1, false, true, true, false, false);

    float techbase = max(stream, jack);
    float techorig = tech;
    tech = (tech / techbase)*tech;
    CalcClamp(tech, techorig*0.85f, techorig);

    float stam;
    if (stream > tech || js > tech || hs > tech)
        if (stream > js && stream > hs)
            stam = Chisel(stream - 0.1f, 2.56f, 1, true, false, true, false, false);
        else if (js > hs)
            stam = Chisel(js - 0.1f, 2.56f, 1, true, false, true, true, false);
        else
            stam = Chisel(hs - 0.1f, 2.56f, 1, true, false, true, false, true);
    else
        stam = Chisel(tech - 0.1f, 2.56f, 1, true, false, false, false, false);

    output.emplace_back(0.f); //temp
    output.emplace_back(downscale_low_accuracy_scores(stream, Scoregoal));

    js = normalizer(js, stream, 7.25f, 0.25f);
    output.emplace_back(downscale_low_accuracy_scores(js, Scoregoal));
    hs = normalizer(hs, stream, 6.5f, 0.3f);
    hs = normalizer(hs, js, 11.5f, 0.15f);
    output.emplace_back(downscale_low_accuracy_scores(hs, Scoregoal));

    float stambase = max(max(stream, tech* 0.96f), max(js, hs));
    if (stambase == stream)
        stambase *= 0.975f;

    stam = normalizer(stam, stambase, 7.75f, 0.2f);
    output.emplace_back(downscale_low_accuracy_scores(stam, Scoregoal));

    output.emplace_back(downscale_low_accuracy_scores(jack, Scoregoal));
    float chordjack = jack * 0.75f;
    output.emplace_back(downscale_low_accuracy_scores(chordjack, Scoregoal));
    float technorm = max(max(stream, js), hs);
    tech = normalizer(tech, technorm, 8.f, .15f) * techscaler;
    output.emplace_back(downscale_low_accuracy_scores(tech, Scoregoal));

    float definitelycj = qprop + hprop + jprop + 0.2f;
    CalcClamp(definitelycj, 0.5f, 1.f);

    // chordjack
    float cj = output[3];

    output[1] *= allhandsdownscaler * manyjumpsdownscaler * lotquaddownscaler;
    output[2] *= nojumpsdownscaler * allhandsdownscaler * lotquaddownscaler;
    output[3] *= nohandsdownscaler * allhandsdownscaler * 1.015f * manyjumpsdownscaler * lotquaddownscaler;
    output[4] *= shortstamdownscaler * 0.985f * lotquaddownscaler;

    cj = normalizer(cj, output[3], 5.5f, 0.3f) * definitelycj * 1.025f;

    if (cj > output[5])
        output[6] = cj;
    else
        output[6] *= 0.9f;

    output[7] *= allhandsdownscaler * manyjumpsdownscaler * lotquaddownscaler * 1.01f;

    float stamclamp = max(max(output[1], output[5]), max(output[2], output[3]));
    CalcClamp(output[4], 1.f, stamclamp * 1.1f);

    dumbvalue = (dumbvalue / static_cast<float>(dumbcounter));
    float stupidvalue = 1.f - (dumbvalue - 2.55f);
    CalcClamp(stupidvalue, 0.85f, 1.f);
    output[7] *= stupidvalue;

    if (stupidvalue <= 0.95f) {
        output[5] *= 1.f + (1.f - sqrt(stupidvalue));
    }

    float skadoot = max(output[3], output[2]);
    if (output[1] < skadoot)
        output[1] -= sqrt(skadoot - output[1]);

    float overall = AggregateScores(output, 0.f, 10.24f, 1);
    output[0] = downscale_low_accuracy_scores(overall, Scoregoal);

    float aDvg = mean(output) * 1.2f;
    for (size_t i = 0; i < output.size(); i++) {
        if (i == 1 || i == 2 || i == 7) {
            CalcClamp(output[i], 0.f, aDvg * 1.0416f);
            output[i] *= grindscaler * grindscaler2;
        } else {
            CalcClamp(output[i], 0.f, aDvg);
            output[i] *= grindscaler * grindscaler2;
        }
        output[i] = downscale_low_accuracy_scores(output[i], Scoregoal);
    }

    output[2] *= jumpthrill;
    output[3] *= jumpthrill;
    output[4] *= sqrt(jumpthrill) * 0.996f;
    output[7] *= sqrt(jumpthrill);

    float highest = 0.f;
    for (auto v : output) {
        if (v > highest)
            highest = v;
    }
    output[0] = AggregateScores(output, 0.f, 10.24f, 1);;

    float dating = 0.5f + (highest / 100.f);
    CalcClamp(dating, 0.f, 0.9f);

    DifficultyRating difficulty = DifficultyRating {output[0],
                                                    output[1],
                                                    output[2],
                                                    output[3],
                                                    output[4],
                                                    output[5],
                                                    output[6],
                                                    output[7]
    };

    if (Scoregoal < dating) {
        difficulty = DifficultyRating {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
    }

    difficulty.jack *= 1.0075f;

    float hsnottech = difficulty.technical - difficulty.handstream;
    float jsnottech = difficulty.technical - difficulty.jumpstream;

    if (highest == difficulty.technical) {
        hsnottech = 4.5f - hsnottech;
        CalcClamp(hsnottech, 0.f, 4.5f);
        difficulty.technical -= hsnottech;

        jsnottech = 4.5f - jsnottech;
        CalcClamp(jsnottech, 0.f, 4.5f);
        difficulty.technical -= jsnottech;
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
        CalcClamp(mod, 1.f, ceil*sqrt(floor));
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
    CalcClamp(output, 0.f, 10000.f);
    return output;
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
            timestamp = (mats1 + mats2 + mats3) / 3.f;

            CalcClamp(timestamp, 25.f, mats3*1.4f);
            float tmp = 1 / timestamp * 2800.f;
            CalcClamp(tmp, 0.f, 50.f);
            output.emplace_back(tmp);
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

    left = new Hand;
    left->InitHand(left_fingers[0], left_fingers[1]);
    left->ohjumpscale = OHJumpDownscaler(NoteInfo, 0, 1);
    left->anchorscale = Anchorscaler(NoteInfo, 0, 1);
    left->rollscale = RollDownscaler(left_fingers[0], left_fingers[1]);
    left->hsscale = HSDownscaler(NoteInfo);
    left->jumpscale = JumpDownscaler(NoteInfo);

    right = new Hand;
    right->InitHand(right_fingers[0], right_fingers[1]);
    right->ohjumpscale = OHJumpDownscaler(NoteInfo, 2, 3);
    right->anchorscale = Anchorscaler(NoteInfo, 2, 3);
    right->rollscale = RollDownscaler(right_fingers[0], right_fingers[1]);
    right->hsscale = left->hsscale;
    right->jumpscale = left->jumpscale;

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
    float Timestamp;
    vector<int> itvnervtmp;
    vector<vector<int>> itvnerv(numitv);

    int left = 1;
    int down = 1 << 1;
    int up = 1 << 2;
    int right = 1 << 3;

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
            Timestamp = 1000 * (scaledtime - last);
            last = scaledtime;
            CalcClamp(Timestamp, 40.f, 5000.f);
            CurrentInterval.emplace_back(Timestamp);
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
    for (size_t i = 0; i < left->v_itvpoints.size(); i++)
        MaxPoints += static_cast<int>(left->v_itvpoints[i] + right->v_itvpoints[i]);
}

float Calc::Chisel(float pskill, float res, int iter, bool stam, bool jack, bool nps, bool js, bool hs) {
    float gotpoints;
    do {
        if (pskill > 100.f)
            return pskill;
        pskill += res;
        if (jack) {
            gotpoints = MaxPoints - JackLoss(j0, pskill) - JackLoss(j1, pskill) - JackLoss(j2, pskill) - JackLoss(j3, pskill);
        }
        else
            gotpoints = left->CalcInternal(pskill, stam, nps, js, hs) + right->CalcInternal(pskill, stam, nps, js, hs);

    } while (gotpoints / MaxPoints < Scoregoal);
    if (iter == 7)
        return pskill;
    return Chisel(pskill - res, res / 2.f, iter + 1, stam, jack, nps, js, hs);
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
        CalcClamp(mod, floor, ceil);
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
vector<float> Calc::OHJumpDownscaler(const vector<NoteInfo>& NoteInfo, int t1, int t2) {
    vector<float> output(nervIntervals.size());
    int firstNote = 1 << t1;
    int secondNote = 1 << t2;

    for (size_t i = 0; i < nervIntervals.size(); i++) {
        if (nervIntervals[i].empty())
            output[i] = 1.f;
        else {
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
    }

    if (SmoothPatterns)
        PatternSmooth(output);
    return output;
}

// pattern modifiers
vector<float> Calc::Anchorscaler(const vector<NoteInfo>& NoteInfo, int t1, int t2) {
    vector<float> output(nervIntervals.size());
    int firstNote = 1 << t1;
    int secondNote = 1 << t2;

    for (size_t i = 0; i < nervIntervals.size(); i++) {
        if (nervIntervals[i].empty())
            output[i] = 1.f;
        else {
            int lcol = 0;
            int rcol = 0;
            for (int row : nervIntervals[i]) {
                if (NoteInfo[row].notes & firstNote)
                    ++lcol;
                if (NoteInfo[row].notes & secondNote)
                    ++rcol;
            }
            bool anyzero = lcol == 0 || rcol == 0;
            output[i] = anyzero ? 1.f : sqrt(1 - (static_cast<float>(min(lcol, rcol)) / static_cast<float>(max(lcol, rcol)) / 4.45f));

            float stupidthing = (static_cast<float>(max(lcol, rcol)) + 2.f) / (static_cast<float>(min(lcol, rcol)) + 1.f);
            dumbvalue += stupidthing;
            ++dumbcounter;

            CalcClamp(output[i], 0.8f, 1.05f);

            if (logpatterns)
                cout << "an " << output[i] << endl;
        }
    }

    if (SmoothPatterns)
        PatternSmooth(output);
    return output;
}


vector<float> Calc::HSDownscaler(const vector<NoteInfo>& NoteInfo) {
    vector<float> output(nervIntervals.size());
    int left = 1;
    int down = 1 << 1;
    int up = 1 << 2;
    int right = 1 << 3;

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
    int down = 1 << 1;
    int up = 1 << 2;
    int right = 1 << 3;

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
            CalcClamp(output[i], 0.f, 1.075f);

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

    left->ohjumpscale.swap(l1);
    left->anchorscale.swap(l2);
    left->rollscale.swap(l3);
    left->hsscale.swap(l4);
    left->jumpscale.swap(l5);

    left->ohjumpscale.shrink_to_fit();
    left->anchorscale.shrink_to_fit();
    left->rollscale.shrink_to_fit();
    left->hsscale.shrink_to_fit();
    left->jumpscale.shrink_to_fit();

    vector<float> r1;
    vector<float> r2;
    vector<float> r3;
    vector<float> r4;
    vector<float> r5;

    right->ohjumpscale.swap(l1);
    right->anchorscale.swap(l2);
    right->rollscale.swap(l3);
    right->hsscale.swap(l4);
    right->jumpscale.swap(l5);

    right->ohjumpscale.shrink_to_fit();
    right->anchorscale.shrink_to_fit();
    right->rollscale.shrink_to_fit();
    right->hsscale.shrink_to_fit();
    right->jumpscale.shrink_to_fit();


    SAFE_DELETE(left);
    SAFE_DELETE(right);
}

// Function to generate SSR rating
DifficultyRating MinaSDCalc(const vector<NoteInfo>& NoteInfo, float musicrate, float goal) {
    unique_ptr<Calc> doot = make_unique<Calc>();
    doot->MusicRate = musicrate;
    CalcClamp(goal, 0.f, 0.965f);	// cap SSR at 96% so things don't get out of hand
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
    {
        for (int i = lower_rate; i < upper_rate; i++)
        {
            allrates.emplace_back(MinaSDCalc(NoteInfo,i / 10.f, 0.93f));
        }
    }
    else
    {
        DifficultyRating output{0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};

        for (int i = lower_rate; i < upper_rate; i++)
        {
            allrates.emplace_back(output);
        }
    }
    return allrates;
}

int GetCalcVersion()
{
    return 263;
}