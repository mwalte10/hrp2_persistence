# ANALYSIS COMPONENT: OUTPUTS (function library: LaTeX tables)
#
# Shared LaTeX table format for every table results_numbers.Rmd writes.
#
# One house style so the generated tables match each other and the
# hand-written appendix tables: caption above, \small text, 4pt column gaps,
# no outer column padding (@{...@{}}, so the rule spans exactly the text),
# booktabs rules, and "median (lo--hi)" on one line. Widths are checked
# against the manuscript's A4 page (12pt, 2.5 cm margins: \textwidth = 160 mm);
# keep new columns short or move detail into the caption.

TAB_FONT   <- "\\small"
TAB_COLSEP <- "4pt"

# Decimal places by quantity, so a value reads the same in every table and in
# the text: half-lives and EIR to 1 dp, sigma_P to 2 dp.
TAB_DP <- c(halflife = 1, eir = 1, sigma = 2)

# "median (lo--hi)" from c(median, lo, hi)
tab_ci <- function(q, d) sprintf("%s (%s--%s)", formatC(q[[1]], format = "f", digits = d, big.mark = ","),
                                 formatC(q[[2]], format = "f", digits = d, big.mark = ","),
                                 formatC(q[[3]], format = "f", digits = d, big.mark = ","))

# Italic group heading spanning all n columns (Prior / Posterior, Population / Site)
tab_group <- function(label, n) sprintf("\\multicolumn{%d}{@{}l}{\\textit{%s}} \\\\", n, label)

# Bold group heading with an optional value in the remaining columns (country rows)
tab_bold_group <- function(label, n, span = n, rest = character(0)) {
  sprintf("\\multicolumn{%d}{@{}l}{%s}%s \\\\", span, label,
          if (length(rest)) paste0(" & ", paste(rest, collapse = " & ")) else "")
}

# A table row: label then cells, indented with \quad when it sits under a group
tab_row <- function(label, cells, indent = TRUE) {
  sprintf("%s%s & %s \\\\", if (indent) "\\quad " else "", label, paste(cells, collapse = " & "))
}

# Floating table (fits on one page)
tex_float_table <- function(colspec, head, body, caption, label) {
  c("\\begin{table}[htbp]", "\\centering", TAB_FONT,
    sprintf("\\setlength{\\tabcolsep}{%s}", TAB_COLSEP),
    sprintf("\\caption{%s}", caption), sprintf("\\label{%s}", label),
    sprintf("\\begin{tabular}{@{}%s@{}}", colspec),
    "\\toprule", head, "\\midrule", body, "\\bottomrule",
    "\\end{tabular}", "\\end{table}")
}

# Multi-page table. The \begingroup keeps the font size and column gap local;
# the header repeats on every page under a "(continued)" caption. A longtable
# is not a float, so it inherits the manuscript's \doublespacing; \linespread{1}
# sets single spacing to match the floating tables.
tex_long_table <- function(colspec, head, body, caption, label) {
  c("\\begingroup", "\\linespread{1}", TAB_FONT,
    sprintf("\\setlength{\\tabcolsep}{%s}", TAB_COLSEP),
    "\\setlength{\\LTcapwidth}{\\textwidth}",
    sprintf("\\begin{longtable}{@{}%s@{}}", colspec),
    sprintf("\\caption{%s}", caption), sprintf("\\label{%s} \\\\", label),
    "\\toprule", head, "\\midrule", "\\endfirsthead",
    "\\caption[]{(continued)} \\\\",
    "\\toprule", head, "\\midrule", "\\endhead",
    "\\bottomrule", "\\endfoot",
    body, "\\end{longtable}", "\\endgroup")
}
