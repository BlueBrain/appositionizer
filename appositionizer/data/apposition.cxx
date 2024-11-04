#include "data/apposition.h"  // Must be on the first substantial line.

#include "data/section.h"
#include "defines.h"

#include <cassert>
#include <iostream>
#include <math.h>

using nlohmann::json;

namespace apposition {

Apposition::Apposition(const Segment& pre, const Segment& post, const std::array<unsigned char, 4>& extends)
    : _pre_key(pre.id())
    , _post_seg(post)
    , _endApposition(extends) {
    bool is_soma = pre.sectionType() == static_cast<int>(cell::SectionType::soma) or
                   post.sectionType() == static_cast<int>(cell::SectionType::soma);
    _angleDotSign = (not is_soma) and (pre.axis() * post.axis() > 0);
}

bool Apposition::operator<(const Apposition& t) const {
    auto lhs = std::tuple{preCell(),
                          postCell(),
                          preSection(),
                          preSegment(),
                          preStart(),
                          postSection(),
                          postSegment(),
                          postStart()};

    auto rhs = std::tuple{t.preCell(),
                          t.postCell(),
                          t.preSection(),
                          t.preSegment(),
                          t.preStart(),
                          t.postSection(),
                          t.postSegment(),
                          t.postStart()};

    return lhs < rhs;
}


bool Region::Connection::assignPostFrom(const Region::AppositionIt& begin, const Region::AppositionIt& end) {
    for (auto itr = begin; itr != end; ++itr) {
        const auto& Apposition = *itr;
        auto s1 = Apposition.preSegment();
        auto la = (int) Apposition.preStart();
        auto ha = (int) Apposition.preEnd();
        auto lb = (int) Apposition.postStart();
        auto hb = (int) Apposition.postEnd();
        auto ta = pre.percentage;

        // if pre region exist in a already existing pre segment (the inital segment array) we
        // calculate
        // depending from the segments direction (dot product) the correct location of the post
        // segment
        // if exist in a gap we flag the post segment id to -1
        if (pre.segment == s1 && ta >= la && ta <= ha) {
            post.segment = Apposition.post();

            if (Apposition.getSegAngleDotSign())  // positive dot product
                post.percentage = ha == la
                                      ? hb
                                      : lb + (unsigned char) int((float(ta - la) / float(ha - la)) *
                                                                 float(hb - lb));
            else  // negative or zero dot product
                post.percentage = ha == la ? lb
                                           : lb + (unsigned char) int(
                                                      (1.0 - float(ta - la) / float(ha - la)) *
                                                      float(hb - lb));

            return true;
        }
    }
    return false;
}

bool Region::Connection::assignPostFrom(const Gaps& gaps) {
    for (const auto& gap: gaps) {
        auto s1pre = gap.first->preSegment();
        auto s1prePortion = (int) gap.first->preEnd();

        auto s2pre = gap.second->preSegment();
        auto s2prePortion = (int) gap.second->preStart();

        auto s1postPortion = (int) gap.first->postStart();
        if ((int) gap.first->getSegAngleDotSign() == true)  // positive dot product
            s1postPortion = (int) gap.first->postEnd();

        // find in which gap the pre Apposition exist and change the pre and post ids to the
        // first point of the Apposition
        if (s1pre == s2pre && s1pre == pre.segment) {
            if (pre.percentage >= s1prePortion && pre.percentage <= s2prePortion) {
                pre.segment = s1pre;
                pre.percentage = s1prePortion;
                post.segment = gap.first->post();
                post.percentage = s1postPortion;
                return true;
            }
        } else if (s1pre < s2pre) {
            bool condition1 = pre.segment >= s1pre && pre.segment <= s2pre;
            bool condition2 = (pre.segment == s1pre && pre.percentage >= s1prePortion) ||
                              (pre.segment == s2pre && pre.percentage <= s2prePortion);
            bool condition3 = (pre.segment > s1pre && pre.segment < s2pre);

            if ((condition1 && condition2) || condition3) {
                pre.segment = s1pre;
                pre.percentage = s1prePortion;
                post.segment = gap.first->post();
                post.percentage = s1postPortion;
                return true;
            }
        }
    }
    return false;
}

void to_json(json& j, const Apposition& t) {
    j = json{{"key1", t.preKey()},
             {"key2", t.postKey()},
             {"ends", {t.preStart(), t.preEnd(), t.postStart(), t.postEnd()}},
             {"dot", t.getSegAngleDotSign()},
             {"segment2_length", t.getSeg2Length()}};
}

void from_json(const json& j, Apposition& t) {
    auto k1 = j["key1"].get<SegmentId>();
    auto k2 = j["key2"].get<SegmentId>();
    auto l = j["segment2_length"].get<t_coordinate>();
    bool angle = j["dot"];

    std::array<unsigned char, 4> extends(j["ends"]);
    static const std::array<t_coordinate, 3> zero{0, 0, 0};
    std::array<t_coordinate, 3> pos{0, 0, l};
    std::array<t_coordinate, 3> neg{0, 0, -l};

    Segment pre(k1, zero, pos);
    Segment post(k2, zero, angle ? pos : neg);

    t = Apposition(pre, post, extends);
}

void to_json(json& j, const Region& r) {
    j = json{{"index", r.getStartAppositionIndex()},
             {"count", r.getSegmentsCount()},
             {"preStart", r.getStartPrePortion()},
             {"preEnd", r.getEndPrePortion()},
             {"preCell", r.getCell1()},
             {"preSection", r.getSection1()},
             {"postCell", r.getCell2()},
             {"segStart", r.getSeg1Begin()},
             {"segEnd", r.getSeg1End()},
             {"connections", r.connections}};
}

void from_json(const json& j, Region& r) {
    cell::SegmentId k1{}, k2{};
    k1.cell(j["preCell"]);
    k1.section(j["preSection"]);
    k1.segment(j["segStart"]);

    k2.cell(j["postCell"]);

    cell::Segment pre(k1);
    cell::Segment post(k2);

    std::array<unsigned char, 4> extends1{j["preStart"], 255, 0, 0};
    std::array<unsigned char, 4> extends2{0, j["preEnd"], 0, 0};

    Apposition t1(pre, post, extends1);

    k1.segment(j["segEnd"]);
    pre = Segment(k1);

    Apposition t2(pre, post, extends2);

    r = Region(t1, t2, j["index"], j["count"]);
    r.connections.clear();
    for (const auto& cj: j["connections"]) {
        r.connections.push_back(cj);
    }
}

void to_json(json& j, const Region::Connection& c) {
    j = json{{"preSegment", c.pre.segment},
             {"prePercentage", c.pre.percentage},
             {"postSection", c.post.segment.section()},
             {"postSegment", c.post.segment.segment()},
             {"postPercentage", c.post.percentage}};
}

void from_json(const json& j, Region::Connection& c) {
    c.pre.segment = j["preSegment"];
    c.pre.percentage = j["prePercentage"];

    cell::SegmentId key(0, j["postSection"], j["postSegment"], 0);

    c.post.segment = cell::Segment(key);
    c.post.percentage = j["postPercentage"];
}

}  // namespace apposition
