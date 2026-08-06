#include <Rcpp.h>
#include <cmath>
#include <vector>

// [[Rcpp::export]]
Rcpp::List rc_corda2_scan_flux_cpp(
    Rcpp::NumericVector flux,
    Rcpp::IntegerVector class_code,
    Rcpp::IntegerVector track_code,
    double threshold
) {
    const R_xlen_t n = flux.size();
    if (class_code.size() != n) {
        Rcpp::stop("CORDA2 flux and class vectors must have equal length.");
    }
    if (!R_FINITE(threshold) || threshold <= 0) {
        Rcpp::stop("CORDA2 flux threshold must be positive and finite.");
    }

    bool tracked[5] = {false, false, false, false, false};
    for (R_xlen_t i = 0; i < track_code.size(); ++i) {
        const int code = track_code[i];
        if (code < 1 || code > 4) {
            Rcpp::stop("CORDA2 tracked class codes must be in 1..4.");
        }
        tracked[code] = true;
    }

    std::vector<int> active;
    std::vector<int> used;
    active.reserve(static_cast<std::size_t>(n / 8 + 1));
    used.reserve(static_cast<std::size_t>(n / 16 + 1));

    for (R_xlen_t i = 0; i < n; ++i) {
        const double value = flux[i];
        if (!R_FINITE(value) || value <= threshold) continue;

        active.push_back(static_cast<int>(i + 1));
        const int code = class_code[i];
        if (code < 1 || code > 4) {
            Rcpp::stop("CORDA2 directional class codes must be in 1..4.");
        }
        if (tracked[code]) {
            used.push_back(static_cast<int>(i + 1));
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("active") = active,
        Rcpp::Named("used") = used
    );
}
