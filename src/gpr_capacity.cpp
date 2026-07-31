#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace {

double and_capacity(
    const Rcpp::NumericMatrix& score,
    int pool,
    int begin,
    int end,
    const Rcpp::IntegerVector& gene_index,
    const std::string& method
) {
    if (begin >= end) return NA_REAL;
    std::vector<double> values;
    values.reserve(static_cast<std::size_t>(end - begin));
    for (int position = begin; position < end; ++position) {
        const int gene = gene_index[position] - 1;
        if (gene < 0 || gene >= score.nrow()) return NA_REAL;
        const double value = score(gene, pool);
        if (!R_FINITE(value)) return NA_REAL;
        values.push_back(value);
    }
    if (method == "min") {
        return *std::min_element(values.begin(), values.end());
    }
    if (method == "mean") {
        long double total = 0.0L;
        for (double value : values) total += value;
        return static_cast<double>(total / values.size());
    }
    std::sort(values.begin(), values.end());
    const std::size_t size = values.size();
    if (size % 2U == 1U) return values[size / 2U];
    return (values[size / 2U - 1U] + values[size / 2U]) / 2.0;
}

double or_capacity(
    const std::vector<double>& values,
    const std::string& method
) {
    if (values.empty()) return NA_REAL;
    if (method == "max") {
        return *std::max_element(values.begin(), values.end());
    }
    if (method == "prob_or") {
        long double product = 1.0L;
        for (double value : values) {
            const double bounded = std::min(1.0, std::max(0.0, value));
            product *= 1.0L - bounded;
        }
        return static_cast<double>(1.0L - product);
    }
    long double total = 0.0L;
    for (double value : values) total += value;
    if (method == "sum_sqrtK") {
        return static_cast<double>(total / std::sqrt(values.size()));
    }
    return static_cast<double>(total);
}

}  // namespace

// [[Rcpp::export]]
Rcpp::NumericMatrix rc_gpr_capacity_cpp(
    Rcpp::NumericMatrix score,
    Rcpp::IntegerVector reaction_group_offset,
    Rcpp::IntegerVector group_gene_offset,
    Rcpp::IntegerVector gene_index,
    std::string and_method,
    std::string or_method
) {
    const int reactions = reaction_group_offset.size() - 1;
    const int groups = group_gene_offset.size() - 1;
    if (reactions < 0 || groups < 0 ||
        reaction_group_offset[0] != 0 || group_gene_offset[0] != 0 ||
        reaction_group_offset[reactions] != groups ||
        group_gene_offset[groups] != gene_index.size()) {
        Rcpp::stop("Compiled GPR offsets are not internally consistent.");
    }
    if (and_method != "min" && and_method != "median" &&
        and_method != "mean") {
        Rcpp::stop("Unsupported GPR AND method.");
    }
    if (or_method != "max" && or_method != "sum_sqrtK" &&
        or_method != "prob_or" && or_method != "sum") {
        Rcpp::stop("Unsupported GPR OR method.");
    }

    Rcpp::NumericMatrix output(reactions, score.ncol());
    std::fill(output.begin(), output.end(), NA_REAL);
    for (int reaction = 0; reaction < reactions; ++reaction) {
        const int group_begin = reaction_group_offset[reaction];
        const int group_end = reaction_group_offset[reaction + 1];
        for (int pool = 0; pool < score.ncol(); ++pool) {
            std::vector<double> available;
            available.reserve(static_cast<std::size_t>(group_end - group_begin));
            for (int group = group_begin; group < group_end; ++group) {
                const double value = and_capacity(
                    score, pool,
                    group_gene_offset[group], group_gene_offset[group + 1],
                    gene_index, and_method
                );
                if (R_FINITE(value)) available.push_back(value);
            }
            output(reaction, pool) = or_capacity(available, or_method);
        }
        if ((reaction & 255) == 0) Rcpp::checkUserInterrupt();
    }
    return output;
}
