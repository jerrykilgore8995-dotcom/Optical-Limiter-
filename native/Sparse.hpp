#pragma once
#include "StampInline.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <bit>
#include <limits>
#include <ostream>
#include <string_view>
#include <vector>

namespace optical {
#if __has_include("GeneratedSolve.hpp")
#include "GeneratedSolve.hpp"
#endif

// Choose and factor the symbolic pattern at prepareToPlay. The audio thread
// reuses a fixed elimination schedule with no searches, allocation or fill-in.
class Sparse {
public:
    void resize(int n) {
        size = n;
        lookup.assign(n * n, -1);
        dense.assign(n * n, 0);
        mutablePattern.assign(n * n, false);
        dynamicPattern.assign(n * n, false);
        recordDynamic = false;
        rhs.assign(n, 0);
        x.assign(n, 0);
        ready = false;
    }

    OPTICAL_STAMP_INLINE void add(int row, int column, double value, bool variable=false) {
        if (row < 0 || column < 0) return;
        const int offset = row * size + column;
        if (!ready) {
            lookup[offset] = 0;
            dense[offset] += value;
            if (variable) mutablePattern[offset]=true;
            if (recordDynamic) dynamicPattern[offset]=true;
        } else {
            const int slot = lookup[offset];
            if (slot >= 0) entries[slot] += value;
        }
    }

    bool prepare() {
        rows.clear();
        updates.clear();
        activePivots.clear();
        std::vector<bool> liveRows(size, true), liveColumns(size, true);
        original.clear();
        for (int row=0; row<size; ++row)
            for (int col=0; col<size; ++col)
                if (has(row,col)) original.push_back({row,col});
        pivoted.assign(size*size,0.);
        pivotRhs.resize(size);
        pivotOrder.resize(size);
        used.resize(size);
        nonzero.resize((size+63)/64*size);
        for (int step = 0; step < size; ++step) {
            std::vector<int> rowCount(size, 0), colCount(size, 0);
            std::vector<double> colMax(size, 0), rowMax(size, 0);
            for (int row = 0; row < size; ++row) if (liveRows[row]) {
                for (int col = 0; col < size; ++col) if (liveColumns[col] && has(row,col)) {
                    ++rowCount[row];
                    ++colCount[col];
                    colMax[col] = std::max(colMax[col], std::abs(v(row,col)));
                    rowMax[row] = std::max(rowMax[row], std::abs(v(row,col)));
                }
            }
            int pr = -1, pc = -1, best = std::numeric_limits<int>::max();
            double tie = -1.0;
            for (int row = 0; row < size; ++row) if (liveRows[row]) {
                for (int col = 0; col < size; ++col) if (liveColumns[col] && has(row,col)) {
                    const double pivot = std::abs(v(row,col));
                    if (pivot < 1e-17 || pivot < 0.01 * colMax[col]) continue;
                    // Prefer self-admittances to off-diagonal transistor
                    // gains which can disappear as the device cuts off.
                    const int cost = (rowCount[row]-1) * (colCount[col]-1) + (row!=col?1:0)
                        + (mutablePattern[row*size+col]?24:0);
                    const double scaled = pivot / (rowMax[row] + 1e-30);
                    if (cost < best || (cost == best && scaled > tie)) {
                        best = cost;
                        tie = scaled;
                        pr = row;
                        pc = col;
                    }
                }
            }
            if (pr < 0) return false;
            liveRows[pr] = liveColumns[pc] = false;
            Row s{};
            s.row = pr;
            s.column = pc;
            for (int col = 0; col < size; ++col)
                if (liveColumns[col] && has(pr,col)) s.upper.push_back({col, 0});
            for (int row = 0; row < size; ++row) {
                if (!liveRows[row] || !has(row,pc)) continue;
                s.lower.push_back({row, 0});
                const double factor = v(row,pc)/v(pr,pc);
                const bool changes = mutablePattern[row*size+pc] || mutablePattern[pr*size+pc];
                for (const auto& upper : s.upper) {
                    lookup[row*size + upper.index] = 0;
                    if (changes || mutablePattern[pr*size+upper.index])
                        mutablePattern[row*size+upper.index]=true;
                    v(row,upper.index) -= factor * v(pr,upper.index);
                }
            }
            rows.push_back(std::move(s));
        }
        int count = 0;
        for (auto& slot : lookup) if (slot >= 0) slot = count++;
        entries.assign(count+1, 0);
        plan.clear();
        uppers.clear();
        lowers.clear();
        for (auto& row : rows) {
            const int upperBegin=static_cast<int>(uppers.size());
            const int lowerBegin=static_cast<int>(lowers.size());
            for (auto& upper : row.upper) {
                upper.slot = at(row.row,upper.index);
                uppers.push_back(upper);
            }
            for (auto& lower : row.lower) {
                lower.slot = at(lower.index,row.column);
                lowers.push_back(lower);
                for (const auto& upper : row.upper)
                    updates.push_back(at(lower.index,upper.index));
            }
            plan.push_back({row.row,row.column,at(row.row,row.column),upperBegin,
                            static_cast<int>(uppers.size()),lowerBegin,static_cast<int>(lowers.size())});
        }
        rows.clear();
        dense.clear();
        ready = true;
        // Default to fully dynamic factorization until a caller supplies a fixed linear part.
        std::vector<int> all(entries.size());
        for (int i=0;i<static_cast<int>(all.size());++i) all[i]=i;
        cacheLinearPart(all);
        return true;
    }

    void cacheLinearPart(const std::vector<int>& mutableSlots) {
        std::vector<bool> dynamic(entries.size(),false);
        for (const int slot : mutableSlots) dynamic[slot]=true;
        runtime.clear();
        activePivots.clear();
        activeUpdates.clear();
        int nextUpdate=0;
        for (const auto& s : plan) {
            const bool pivotChanges=dynamic[s.pivot];
            const double reciprocal=pivotChanges?0:1/entries[s.pivot];
            const int begin=static_cast<int>(runtime.size());
            for (int j=s.lowerBegin;j<s.lowerEnd;++j) {
                const auto& row=lowers[j];
                const bool factorChanges=pivotChanges||dynamic[row.slot];
                const double factor=factorChanges?0:entries[row.slot]*reciprocal;
                if (!factorChanges && factor == 0.) {
                    nextUpdate += s.upperEnd - s.upperBegin;
                    continue;
                }
                const int first=static_cast<int>(activeUpdates.size());
                for (int k=s.upperBegin;k<s.upperEnd;++k) {
                    const int target=updates[nextUpdate++];
                    const int source=uppers[k].slot;
                    if (!factorChanges && !dynamic[source]) entries[target]-=factor*entries[source];
                    else {
                        dynamic[target]=true;
                        activeUpdates.push_back({target,source});
                    }
                }
                runtime.push_back({row.index,row.slot,first,static_cast<int>(activeUpdates.size()),
                                   factor,factorChanges});
            }
            activePivots.push_back({s.row,s.pivot,begin,static_cast<int>(runtime.size()),
                                    reciprocal,pivotChanges});
        }
        fixed=entries;
        signature=0xcbf29ce484222325ULL;
        const auto add=[&](std::uint64_t value) {
            signature^=value;
            signature*=0x100000001b3ULL;
        };
        for (const auto& row : activePivots) {
            for (int v : {row.row,row.pivot,row.begin,row.end}) add(v);
            add(row.changes);
            add(std::bit_cast<std::uint64_t>(row.reciprocal));
        }
        for (const auto& row : runtime) {
            for (int v : {row.index,row.slot,row.first,row.last}) add(v);
            add(row.changes);
            add(std::bit_cast<std::uint64_t>(row.factor));
        }
        for (const auto& e : activeUpdates) {add(e.index);add(e.slot);}
        for (const auto& s : plan) {
            for (int v : {s.row,s.column,s.pivot,s.upperBegin,s.upperEnd}) add(v);
        }
        for (const auto& e : uppers) {add(e.index);add(e.slot);}
#if __has_include("GeneratedSolve.hpp")
        generated=findGeneratedSolve(signature);
#endif
    }

    void startIteration(const std::vector<double>& initial) {
        std::copy(fixed.begin(),fixed.end(),entries.begin());
        std::copy(initial.begin(),initial.end(),rhs.begin());
    }

    void saveLinearEntries() { linearEntries = entries; }

    void beginDynamicStamp() {
        recordDynamic = true;
    }

    void cachePreparedFactorization() {
        recordDynamic = false;
        std::vector<int> dynamicSlots;
        for (const auto& position : original)
            if (dynamicPattern[position.index * size + position.slot])
                dynamicSlots.push_back(at(position.index, position.slot));
        restoreLinearEntries(true);
        cacheLinearPart(dynamicSlots);
        std::vector<bool>{}.swap(dynamicPattern);
    }

    void restoreLinearEntries(bool raw = false) {
        const auto& initial = raw ? linearEntries : fixed;
        std::copy(initial.begin(), initial.end(), entries.begin());
    }

    void clear() {
        std::fill(entries.begin(), entries.end(), 0.0);
        std::fill(rhs.begin(), rhs.end(), 0.0);
    }

    void emit(std::ostream& code,std::string_view name) const {
        code<<"inline constexpr std::uint64_t signature_"<<name<<'='<<signature<<"ULL;\n";
        code<<"inline bool solve_"<<name<<"(double* __restrict a,double* __restrict b,double* __restrict x) {\n";
        code<<"    double u=0,v=0;\n";
        for (const auto& s : activePivots) {
            if (s.changes) {
                code<<"    u=a["<<s.pivot<<"]; if (std::abs(u)<1e-17 || !std::isfinite(u)) return false; "
                    <<"u=1/u; a["<<s.pivot<<"]=u;\n";
            }
            for (int j=s.begin;j<s.end;++j) {
                const auto& row=runtime[j];
                if (row.changes) {
                    code<<"    v=a["<<row.slot<<"]*";
                    if (s.changes) code<<"u";
                    else code<<std::hexfloat<<s.reciprocal<<std::defaultfloat;
                    code<<";\n";
                }
                else if (row.factor == 0) continue;
                else code<<"    v="<<std::hexfloat<<row.factor<<std::defaultfloat<<";\n";
                code<<"    b["<<row.index<<"]-=v*b["<<s.row<<"];\n";
                for (int k=row.first;k<row.last;++k)
                    code<<"    a["<<activeUpdates[k].index<<"]-=v*a["<<activeUpdates[k].slot<<"];\n";
            }
        }
        for (auto s=plan.rbegin();s!=plan.rend();++s) {
            code<<"    v=b["<<s->row<<"];\n";
            for (int k=s->upperBegin;k<s->upperEnd;++k)
                code<<"    v-=a["<<uppers[k].slot<<"]*x["<<uppers[k].index<<"];\n";
            const auto& forward=activePivots[static_cast<std::size_t>((s+1).base()-plan.begin())];
            if (forward.changes) code<<"    x["<<s->column<<"]=v*a["<<s->pivot<<"];\n";
            else code<<"    x["<<s->column<<"]=v*"<<std::hexfloat<<forward.reciprocal
                     <<std::defaultfloat<<";\n";

        }
        code<<"    bool finite=true;\n";
        code<<"    for (int i=0;i<"<<size<<";++i) finite &= std::isfinite(x[i]);\n";
        code<<"    return finite;\n}\n";
    }

    bool solve() {
#if __has_include("GeneratedSolve.hpp")
        if (generated) return generated(entries.data(),rhs.data(),x.data());
#endif
        const auto* upper=uppers.data();
        auto* values=entries.data();
        auto* vector=rhs.data();
        for (const auto& s : activePivots) {
            double reciprocal=s.reciprocal;
            if (s.changes) {
                const double pivot=values[s.pivot];
                if (std::abs(pivot)<1e-17 || !std::isfinite(pivot)) return false;
                reciprocal=1.0/pivot;
            }
            values[s.pivot]=reciprocal;
            const double y=vector[s.row];
            for (int j=s.begin;j<s.end;++j) {
                const auto& row=runtime[j];
                const double factor=row.changes?values[row.slot]*reciprocal:row.factor;
                vector[row.index]-=factor*y;
                for (int k=row.first;k<row.last;++k)
                    values[activeUpdates[k].index]-=factor*values[activeUpdates[k].slot];
            }
        }
        for (auto s=plan.rbegin();s!=plan.rend();++s) {
            double y=vector[s->row];
            for (int k=s->upperBegin;k<s->upperEnd;++k)
                y-=values[upper[k].slot]*x[upper[k].index];
            x[s->column]=y*values[s->pivot];
            if (!std::isfinite(x[s->column])) return false;
        }
        return true;
    }

    // A nonlinear cutoff can invalidate a formerly sound static pivot.
    // Refactor the same Jacobian with numeric row pivoting in that case.
    bool solvePivoted() {
        const int words=(size+63)/64;
        std::fill(pivoted.begin(),pivoted.end(),0.);
        std::fill(nonzero.begin(),nonzero.end(),0);
        std::fill(used.begin(),used.end(),false);
        pivotRhs=rhs;
        for (const auto& p : original) {
            const double value=entries[at(p.index,p.slot)];
            if (value!=0) {
                pivoted[p.index*size+p.slot]=value;
                nonzero[p.index*words+p.slot/64]|=std::uint64_t{1}<<(p.slot%64);
            }
        }
        for (int i=0;i<size;++i) {
            const int col=plan[i].column;
            int row=-1;
            double best=0.;
            for (int candidate=0;candidate<size;++candidate) {
                if (used[candidate]) continue;
                const double value=std::abs(pivoted[candidate*size+col]);
                if (value>best) {best=value;row=candidate;}
            }
            if (row<0 || best<1e-18 || !std::isfinite(best)) return false;
            const int preferred=plan[i].row;
            if (!used[preferred] && std::abs(pivoted[preferred*size+col])>best*.1)
                row=preferred;
            used[row]=true;
            pivotOrder[i]=row;
            auto* upper=pivoted.data()+row*size;
            const double main=upper[col];
            for (int target=0;target<size;++target) {
                auto* lower=pivoted.data()+target*size;
                if (used[target] || lower[col]==0.) continue;
                const double factor=lower[col]/main;
                pivotRhs[target]-=factor*pivotRhs[row];
                lower[col]=0;
                auto* mask=nonzero.data()+target*words;
                mask[col/64]&=~(std::uint64_t{1}<<(col%64));
                for (int word=0;word<words;++word) {
                    auto bits=nonzero[row*words+word];
                    for (;bits;bits&=bits-1) {
                        const int index=word*64+std::countr_zero(bits);
                        if (index!=col) lower[index]-=factor*upper[index];
                    }
                    mask[word]|=nonzero[row*words+word];
                }
                mask[col/64]&=~(std::uint64_t{1}<<(col%64));
            }
        }
        for (int i=size-1;i>=0;--i) {
            const int row=pivotOrder[i], col=plan[i].column;
            auto* upper=pivoted.data()+row*size;
            double y=pivotRhs[row];
            for (int word=0;word<words;++word) {
                auto bits=nonzero[row*words+word];
                for (;bits;bits&=bits-1) {
                    const int index=word*64+std::countr_zero(bits);
                    if (index!=col) y-=upper[index]*x[index];
                }
            }
            x[col]=y/upper[col];
            if (!std::isfinite(x[col])) return false;
        }
        return true;
    }

    int slot(int row,int col) const {
        if (row<0 || col<0) return static_cast<int>(entries.size()-1);
        const int i=at(row,col);
        return i<0?static_cast<int>(entries.size()-1):i;
    }
    void addSlot(int slot,double value) { entries[slot]+=value; }

    int operationCount() const { return static_cast<int>(activeUpdates.size()); }
    std::vector<double> rhs, x;

private:
    struct Terminal { int index, slot; };
    struct Row {
        int row, column;
        std::vector<Terminal> upper, lower;
    };
    int at(int row, int col) const { return lookup[row*size+col]; }
    bool has(int row, int col) const { return at(row,col) >= 0; }
    double& v(int row, int col) { return dense[row*size+col]; }

    std::uint64_t signature=0;
#if __has_include("GeneratedSolve.hpp")
    bool (*generated)(double*,double*,double*)=nullptr;
#endif
    int size = 0;
    bool ready = false;
    std::vector<int> lookup, updates;
    std::vector<bool> dynamicPattern;
    bool recordDynamic = false;
    std::vector<Row> rows;
    struct Pivot { int row,column,pivot,upperBegin,upperEnd,lowerBegin,lowerEnd; };
    std::vector<Pivot> plan;
    struct ActiveLower { int index,slot,first,last;double factor;bool changes; };
    struct ActivePivot { int row,pivot,begin,end;double reciprocal;bool changes; };
    std::vector<ActiveLower> runtime;
    std::vector<ActivePivot> activePivots;
    std::vector<Terminal> uppers,lowers,activeUpdates;
    std::vector<double> dense, entries, fixed, linearEntries;
    std::vector<Terminal> original;
    std::vector<double> pivoted,pivotRhs;
    std::vector<int> pivotOrder;
    std::vector<bool> used;
    std::vector<std::uint64_t> nonzero;
    std::vector<bool> mutablePattern;
};
}
