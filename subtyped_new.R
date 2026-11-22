# =====================================================================
# subtyped_new.R
#
# REVISION COPY (2026-11-22): fixes marked FIX 1-8 in the code.
#   1. read_excel(..., guess_max = 10000)        changes results (small)
#   2. start/end years swapped where end < start changes results (Nigeria)
#   3. risk-of-bias score compared as a number   no change now
#   4. CRF07_BC NA no longer wipes out other_CRF  count files only
#   5. (reverted: inv_simpson name kept)
#   6. country-periods with no HIV-1 counts dropped from indices.csv
#   7. RMSD and diversity indices over classified variants only (k = 172);
#      % unclassified added to indices.csv          changes results
#   8. sample size and % recombinant / CRF / URF use the summed variant
#      counts (total_recorded_genotypes_adj), not the raw "Total number
#      genotyped" field                              changes results
# CHECKS (2026-11-22): lines starting "# CHECK" test that the counts add
# up at each step. stopifnot() stops the script if something that should
# never happen does; message() lines print a summary to look over.
#
# Takes the REDCap export of the systematic review and produces:
#   counts_country_subtype_distirbutions.csv  adjusted variant counts by
#                                             country and time period
#   proportions.csv                           the same, as proportions
#   RMSDcountries.csv                         pairwise RMSD between every
#                                             country pair in each period
#   indices.csv                               Shannon / Simpson diversity
#                                             and recombinant proportions
#
# Steps
#   1. Load and clean the review data
#   2. Assign regions
#   3. Coerce counts to numeric, split studies across years
#   4. Collapse the variant columns
#   5. Aggregate by country and time period
#   6. Pairwise RMSD
#   7. Diversity indices
# =====================================================================

# Packages are loaded by 00_libraries.R (see run_all.R).

# =====================================================================
# 1. Load and clean
# =====================================================================

# note that this reads in the full dataset. add the suffix '_published' to run with the restricted data
# FIX 1 guess_max: read_excel guesses each column's type from the first 1000
# rows by default. Columns that are empty in those rows are read as
# logical and any later values are silently set to NA (this lost
# CRF34_01B x5 and CRF73_BG x1). Guessing from all rows avoids that.
HIV_samples <- read_excel("HIVIDDONDPHProject_DATA_LABELS_2024-05-07_1226.xlsx", guess_max = 10000)
#HIV_samples <- read_excel("raw data/HIVIDDONDPHProject_DATA_LABELS_2024-05-07_1226_published.xlsx")

# FIX 3: the export stores every column as text, so the score was being
# compared as a string ("4.0" < 4). That happens to work for the scores
# present (1-9) but would misclassify 10+; compare it as a number.
HIV_samples <- HIV_samples %>%
  mutate(`Summary score for overall risk of study bias` =
           as.numeric(`Summary score for overall risk of study bias`)) %>%
  mutate(`Summary of overall risk of study bias` = case_when(
    `Summary score for overall risk of study bias` < 4 ~ "Low Risk",
    `Summary score for overall risk of study bias` > 6 ~ "High Risk",
    `Summary score for overall risk of study bias` > 3 &
      `Summary score for overall risk of study bias` < 7 ~ "Moderate Risk")) %>%
  mutate(`Site 1: Country` = ifelse(
    `Site 1: Country` == "United States Minor Outlying Islands (the)",
    "United States of America (the)", `Site 1: Country`))

# Drop records whose sampling sites span more than one country: only
# single-country records can be assigned to a country distribution.
to_throw <- HIV_samples %>%
  filter(!is.na(.data[[paste0("Site ", 2, ": Country")]])) %>%
  filter(if_any(all_of(paste0("Site ", 2, ": Country")),
                ~ . != .data[[paste0("Site 1: Country")]]))

for (i in 3:30) {
  temp <- HIV_samples %>%
    filter(!is.na(.data[[paste0("Site ", i, ": Country")]])) %>%
    filter(if_any(all_of(paste0("Site ", i, ": Country")),
                  ~ . != .data[[paste0("Site 1: Country")]]))
  to_throw <- plyr::rbind.fill(to_throw, temp)
}

to_throw <- to_throw %>% distinct()

HIV_samples <- HIV_samples %>% anti_join(to_throw)

# CHECK: report how many multi-country records were removed
message("CHECK multi-country records removed: ", nrow(to_throw),
        "; records left: ", nrow(HIV_samples))
stopifnot(!any(HIV_samples$`Record Number` %in% to_throw$`Record Number`))
stopifnot(!any(duplicated(HIV_samples$`Record Number`)))


# =====================================================================
# 2. Regions
#

# =====================================================================
caribbean <- c("Bahamas", "Barbados", "Belize", "Cuba", "Dominican Republic",
               "Dominican Republic (the)", "Guadeloupe", "Haiti", "Jamaica",
               "Martinique", "Puerto Rico", "Trinidad and Tobago")

latin_america <- c("Argentina", "Bolivia", "Bolivia (Plurinational State of)",
                   "Brazil", "Chile", "Colombia", "Costa Rica", "Ecuador",
                   "El Salvador", "Guatemala", "Guyana", "French Guiana",
                   "Honduras", "Mexico", "Nicaragua", "Panama", "Paraguay",
                   "Peru", "Suriname", "Uruguay", "Venezuela",
                   "Venezuela (Bolivarian Republic of)")

north_america <- c("United States of America", "United States of America (the)",
                   "United States Minor Outlying Islands (the)", "Canada")

# Western and Central Europe
wce <- c("Austria", "Belgium", "Bulgaria", "Croatia", "Cyprus", "Czech Republic",
         "Czechia", "Denmark", "Estonia", "Finland", "France", "Germany",
         "Greece", "Greenland", "Hungary", "Iceland", "Ireland", "Israel",
         "Italy", "Latvia", "Lithuania", "Luxembourg", "Malta", "Netherlands",
         "Netherlands (the)", "Norway", "Poland", "Portugal", "Romania",
         "Serbia", "Slovakia", "Slovenia", "Spain", "Sweden", "Switzerland",
         "Turkey", "United Kingdom of Great Britain & Northern Ireland",
         "United Kingdom of Great Britain and Northern Ireland (the)")

# Eastern Europe & Central Asia
eeca <- c("Albania", "Armenia", "Azerbaijan", "Belarus", "Bosnia and Herzegovina",
          "Georgia", "Kazakstan", "Kazakhstan", "Kyrgyzstan", "Montenegro",
          "Republic of Moldova", "Moldova (the Republic of)", "Russian Federation",
          "Russian Federation (the)", "Tajikistan",
          "The former Yugoslav Republic of Macedonia", "Ukraine", "Uzbekistan")

india_nepal_sl <- c("India", "Nepal", "Sri Lanka")

southeast_asia <- c("Afghanistan", "Bangladesh", "Bhutan", "Brunei Darussalam",
                    "Cambodia", "Indonesia", "Lao People's Democratic Republic",
                    "Lao People's Democratic Republic (the)", "Malaysia",
                    "Maldives", "Myanmar", "Pakistan", "Iran (Islamic Republic of)",
                    "Philippines", "Philippines (the)", "Singapore", "Thailand",
                    "Timor-Leste", "Viet Nam")

east_asia <- c("China", "Democratic People's Republic of Korea", "Hong Kong",
               "Japan", "Mongolia", "Republic of Korea", "Korea (the Republic of)",
               "Taiwan (Province of China)")

oceania <- c("Australia", "Fiji", "New Zealand", "Papua New Guinea")

middle_east_north_africa <- c("Algeria", "Egypt", "Kuwait", "Lebanon", "Libya",
                              "Morocco", "Oman", "Saudi Arabia", "Tunisia",
                              "Yemen", "W. Sahara")

west_africa <- c("Benin", "Burkina Faso", "Cameroon", "Cape Verde", "Cabo Verde",
                 "Cote d'Ivoire", "Côte d'Ivoire", "CÃ´te d'Ivoire", "Gambia",
                 "Gambia (the)", "Ghana", "Guinea", "Guinea-Bissau", "Liberia",
                 "Mali", "Mauritania", "Niger", "Niger (the)", "Nigeria",
                 "Senegal", "Sierra Leone", "Togo")

central_africa <- c("Angola", "Central African Republic", "Central African Rep.",
                    "Central African Republic (the)", "Chad", "Congo", "Congo (the)",
                    "Democratic Republic of the Congo",
                    "Congo (the Democratic Republic of the)", "Dem. Rep. Congo",
                    "Equatorial Guinea", "Eq. Guinea", "Gabon")

southern_africa <- c("Botswana", "Eswatini", "eSwatini", "Lesotho", "Malawi",
                     "Mozambique", "Namibia", "South Africa", "Swaziland",
                     "Zambia", "Zimbabwe")

eth_erit_dji <- c("Ethiopia", "Eritrea", "Djibouti")

east_africa <- c("Burundi", "Kenya", "Madagascar", "Mauritius", "Rwanda",
                 "Somalia", "South Sudan", "S. Sudan", "Sudan", "Sudan (the)",
                 "Uganda", "United Republic of Tanzania",
                 "Tanzania, the United Republic of", "Tanzania", "Somaliland")

HIV_samples <- HIV_samples %>%
  mutate(region = case_when(
    `Site 1: Country` %in% caribbean                ~ "caribbean",
    `Site 1: Country` %in% latin_america            ~ "latin_america",
    `Site 1: Country` %in% north_america            ~ "north_america",
    `Site 1: Country` %in% wce                      ~ "wce",
    `Site 1: Country` %in% eeca                     ~ "eeca",
    `Site 1: Country` %in% india_nepal_sl           ~ "india_nepal_sl",
    `Site 1: Country` %in% southeast_asia           ~ "southeast_asia",
    `Site 1: Country` %in% east_asia                ~ "east_asia",
    `Site 1: Country` %in% oceania                  ~ "oceania",
    `Site 1: Country` %in% middle_east_north_africa ~ "middle_east_north_africa",
    `Site 1: Country` %in% west_africa              ~ "west_africa",
    `Site 1: Country` %in% east_africa              ~ "east_africa",
    `Site 1: Country` %in% eth_erit_dji             ~ "eth_erit_dji",
    `Site 1: Country` %in% central_africa           ~ "central_africa",
    `Site 1: Country` %in% southern_africa          ~ "southern_africa",
    TRUE ~ NA_character_
  ))

# CHECK: every country has a region (a country missing from the lists
# above would be NA here)
no_region <- HIV_samples %>% filter(is.na(region)) %>% distinct(`Site 1: Country`)
message("CHECK countries with no region: ", nrow(no_region),
        if (nrow(no_region) > 0) paste0(" (", paste(no_region$`Site 1: Country`, collapse = "; "), ")") else "")
stopifnot(nrow(no_region) == 0)


# =====================================================================
# 3. Column types and study years
#
# Studies spanning several years are split across those years, each year
# receiving 1 / no_years of the study's variant counts.
# =====================================================================

# Genotyping fragment columns, found by name rather than by position
# (positions 343:352 in this export), so a change in column order cannot
# silently pick the wrong columns.
column_indices <- grep("^Genotyping fragment", names(HIV_samples))
stopifnot(length(column_indices) == 10)
for (i in column_indices) {
  HIV_samples[[i]] <- as.character(HIV_samples[[i]])
}

subtype_cols <- c(grep("^HIV-1 group [MNOP]", names(HIV_samples), value = TRUE),
                  grep("^CRF", names(HIV_samples), value = TRUE),
                  grep("^URF\\d+ number$", names(HIV_samples), value = TRUE),
                  "Number of unspecified CRFs", "Number of undefined URFs",
                  "Unspecified recombinants", "Total number genotyped")

# CHECK (report only): text in a count column becomes NA below and those
# samples are lost, so list any (record 1016 has "C" where the count of
# subtype C should be)
for (col in intersect(subtype_cols, names(HIV_samples))) {
  not_numeric <- !is.na(HIV_samples[[col]]) & is.na(suppressWarnings(as.numeric(HIV_samples[[col]])))
  if (any(not_numeric)) {
    message("CHECK text in count column '", col, "', set to NA: record(s) ",
            paste(HIV_samples$`Record Number`[not_numeric], collapse = ", "), " (value: ",
            paste(HIV_samples[[col]][not_numeric], collapse = ", "), ")")
  }
}

HIV_samples <- HIV_samples %>%
  mutate(across(any_of(subtype_cols), ~ suppressWarnings(as.numeric(.))))

# CHECK: no negative counts in any variant column
min_count <- min(as.matrix(HIV_samples[, intersect(subtype_cols, names(HIV_samples))]), na.rm = TRUE)
message("CHECK smallest count in any variant column: ", min_count)
stopifnot(min_count >= 0)

# CHECK (report only): counts that are not whole numbers (record 389 has
# 3.7 for subtype A)
count_matrix_raw <- as.matrix(HIV_samples[, intersect(subtype_cols, names(HIV_samples))])
not_whole <- which(!is.na(count_matrix_raw) & abs(count_matrix_raw - round(count_matrix_raw)) > 1e-9, arr.ind = TRUE)
message("CHECK counts that are not whole numbers: ", nrow(not_whole),
        if (nrow(not_whole) > 0) paste0(" (", paste0("record ", HIV_samples$`Record Number`[not_whole[, 1]], " ",
                                                     colnames(count_matrix_raw)[not_whole[, 2]], " = ",
                                                     count_matrix_raw[not_whole], collapse = "; "), ")") else "")

# CHECK: keep the record-level total of the variant counts (every
# HIV-1 variant column, NOT the raw Total number genotyped field), to
# test that splitting across years and aggregating does not gain or
# lose samples
variant_cols_raw <- setdiff(intersect(subtype_cols, names(HIV_samples)), "Total number genotyped")
records_total_variants <- sum(HIV_samples[, variant_cols_raw], na.rm = TRUE)

HIV_samples <- HIV_samples %>%
  mutate(`Year study started`     = as.numeric(`Year study started`),
         `Year study ended`       = as.numeric(`Year study ended`),
         `Total number genotyped` = as.numeric(`Total number genotyped`)) %>%
  mutate(`Year study ended` = ifelse(is.na(`Year study ended`),
                                     `Year study started`, `Year study ended`)) %>%
  # FIX 2: swap start and end where they were entered the wrong way round
  # (record 795, Nigeria, 2009-2006); otherwise no_years is negative and
  # the record's samples are subtracted rather than added
  mutate(start_fixed          = pmin(`Year study started`, `Year study ended`),
         end_fixed            = pmax(`Year study started`, `Year study ended`),
         `Year study started` = start_fixed,
         `Year study ended`   = end_fixed) %>%
  dplyr::select(-start_fixed, -end_fixed) %>%
  mutate(no_years           = `Year study ended` - `Year study started` + 1,
         freq_var           = 1 / no_years,
         year               = map2(`Year study started`, `Year study ended`, seq)) %>%
  unnest(cols = year)

# CHECK: every record spans at least one year, and each record's year
# weights (freq_var) add up to 1, so no samples are gained or lost
stopifnot(all(!is.na(HIV_samples$no_years)), all(HIV_samples$no_years >= 1))
weight_check <- HIV_samples %>%
  group_by(`Record Number`) %>%
  summarise(weight = sum(freq_var), n_rows = n(), no_years = first(no_years))
stopifnot(all(abs(weight_check$weight - 1) < 1e-9),
          all(weight_check$n_rows == weight_check$no_years))
message("CHECK years per record: ", min(weight_check$no_years), " to ",
        max(weight_check$no_years), "; ", sum(weight_check$no_years > 1),
        " of ", nrow(weight_check), " records span more than one year")
split_total_variants <- sum(rowSums(HIV_samples[, variant_cols_raw], na.rm = TRUE) * HIV_samples$freq_var)
stopifnot(abs(split_total_variants - records_total_variants) < 1e-6)
message("CHECK total variant counts after splitting across years: ",
        round(split_total_variants, 1), " (records: ", records_total_variants, ")")

HIV_samples <- HIV_samples %>%
  mutate(year_category = case_when(
    year >= 1980 & year <= 1989 ~ "1980-1989",
    year >= 1990 & year <= 1994 ~ "1990-1994",
    year >= 1995 & year <= 1999 ~ "1995-1999",
    year >= 2000 & year <= 2004 ~ "2000-2004",
    year >= 2005 & year <= 2009 ~ "2005-2009",
    year >= 2010 & year <= 2014 ~ "2010-2014",
    year >= 2015 & year <= 2019 ~ "2015-2019",
    year >  2019                ~ "2020-2022",
    TRUE ~ NA_character_
  ))


# =====================================================================
# 4. Collapse the variant columns
#
# Sub-subtypes are summed into their parent (A1..A8 into A, F1/F2 into F).
# HIV-2 is not analysed.
# =====================================================================
# CHECK: every A and F sub-subtype column in the export is summed below
A_parts <- c("HIV-1 group M: A", "HIV-1 group M: A1", "HIV-1 group M: A2",
             "HIV-1 group M: A3", "HIV-1 group M: A4", "HIV-1 group M: A6",
             "HIV-1 group M: A7", "HIV-1 group M: A8")
F_parts <- c("HIV-1 group M: F", "HIV-1 group M: F1", "HIV-1 group M: F2")
stopifnot(setequal(grep("^HIV-1 group M: A\\d*$", names(HIV_samples), value = TRUE), A_parts),
          setequal(grep("^HIV-1 group M: F\\d*$", names(HIV_samples), value = TRUE), F_parts))
# CHECK: report any group M/N/O/P column that is not used below
used_group_cols <- c(A_parts, F_parts, paste0("HIV-1 group M: ", c("B", "C", "D", "G", "H", "J", "K", "L")),
                     "HIV-1 group N", "HIV-1 group O", "HIV-1 group P")
unused_group_cols <- setdiff(grep("^HIV-1 group", names(HIV_samples), value = TRUE), used_group_cols)
message("CHECK HIV-1 group columns not used: ",
        if (length(unused_group_cols) == 0) "none" else paste(unused_group_cols, collapse = "; "))
stopifnot(length(unused_group_cols) == 0)

HIV_samples$A <- rowSums(HIV_samples[, c("HIV-1 group M: A",  "HIV-1 group M: A1",
                                         "HIV-1 group M: A2", "HIV-1 group M: A3",
                                         "HIV-1 group M: A4", "HIV-1 group M: A6",
                                         "HIV-1 group M: A7", "HIV-1 group M: A8")],
                         na.rm = TRUE)
HIV_samples$B   <- HIV_samples$`HIV-1 group M: B`
HIV_samples$C   <- HIV_samples$`HIV-1 group M: C`
HIV_samples$D   <- HIV_samples$`HIV-1 group M: D`
HIV_samples$`F` <- rowSums(HIV_samples[, c("HIV-1 group M: F", "HIV-1 group M: F1",
                                           "HIV-1 group M: F2")], na.rm = TRUE)
HIV_samples$G <- HIV_samples$`HIV-1 group M: G`
HIV_samples$H <- HIV_samples$`HIV-1 group M: H`
HIV_samples$J <- HIV_samples$`HIV-1 group M: J`
HIV_samples$K <- HIV_samples$`HIV-1 group M: K`
HIV_samples$L <- HIV_samples$`HIV-1 group M: L`
HIV_samples$N <- HIV_samples$`HIV-1 group N`
HIV_samples$O <- HIV_samples$`HIV-1 group O`
HIV_samples$P <- HIV_samples$`HIV-1 group P`

# CRFs: all named CRF columns plus the unspecified count
CRF_cols <- grep("^CRF", names(HIV_samples), value = TRUE)
CRF_cols <- c(CRF_cols, "Number of unspecified CRFs")

HIV_samples$total_CRF <- rowSums(HIV_samples[, c(CRF_cols)], na.rm = TRUE)

# CRF01_AE, CRF02_AG and CRF07_BC are shown individually in the figures,
# so other_CRF must exclude all three.
# FIX 4: CRF07_BC was not replace_na'd, so other_CRF was NA in every
# record without CRF07_BC (1,543 of 1,930 records with CRFs) and the
# other_CRF column held 14,995 instead of 31,905 adjusted samples. The
# pie charts recompute other_CRF themselves, so figures were not
# affected; the other_CRF columns in the count files were.
HIV_samples <- HIV_samples %>%
  mutate(CRF01_AE  = replace_na(CRF01_AE, 0),
         CRF02_AG  = replace_na(CRF02_AG, 0),
         CRF07_BC  = replace_na(CRF07_BC, 0),
         other_CRF = total_CRF - (CRF01_AE + CRF02_AG + CRF07_BC))

# URFs: all numbered URF columns plus the undefined count
matching_cols <- grep("^URF\\d+ number$", names(HIV_samples), value = TRUE)
matching_cols <- c(matching_cols, "Number of undefined URFs")

HIV_samples$URFs <- rowSums(HIV_samples[, c(matching_cols)], na.rm = TRUE)

HIV_samples <- HIV_samples %>%
  mutate(unspecified_recombinants = replace_na(`Unspecified recombinants`, 0))

HIV_samples <- HIV_samples %>%
  mutate(total_recombinants = URFs + total_CRF + unspecified_recombinants)

# CHECK: in every row, the collapsed variant columns add up to the sum of
# all the raw HIV-1 variant columns, so no raw column was dropped or
# counted twice when collapsing; other_CRF is never NA or negative
collapsed_total <- HIV_samples$A + replace_na(HIV_samples$B, 0) + replace_na(HIV_samples$C, 0) +
  replace_na(HIV_samples$D, 0) + HIV_samples$`F` + replace_na(HIV_samples$G, 0) +
  replace_na(HIV_samples$H, 0) + replace_na(HIV_samples$J, 0) + replace_na(HIV_samples$K, 0) +
  replace_na(HIV_samples$L, 0) + replace_na(HIV_samples$N, 0) + replace_na(HIV_samples$O, 0) +
  replace_na(HIV_samples$P, 0) + HIV_samples$total_CRF + HIV_samples$URFs + HIV_samples$unspecified_recombinants
stopifnot(all(abs(collapsed_total - rowSums(HIV_samples[, variant_cols_raw], na.rm = TRUE)) < 1e-9),
          all(!is.na(HIV_samples$other_CRF)), all(HIV_samples$other_CRF >= 0))
message("CHECK collapsed columns: ", length(CRF_cols) - 1, " named CRFs + unspecified CRFs, ",
        length(matching_cols) - 1, " URF columns + undefined URFs")

# CHECK (report only): records whose variant counts add up to MORE than
# Total number genotyped (fewer is expected, e.g. HIV-2 or failed
# subtyping; more suggests a data entry error)
record_variants <- HIV_samples %>%
  distinct(`Record Number`, .keep_all = TRUE) %>%
  mutate(variant_sum = A + replace_na(B, 0) + replace_na(C, 0) + replace_na(D, 0) + `F` +
           replace_na(G, 0) + replace_na(H, 0) + replace_na(J, 0) + replace_na(K, 0) +
           replace_na(L, 0) + replace_na(N, 0) + replace_na(O, 0) + replace_na(P, 0) +
           URFs + total_CRF + unspecified_recombinants) %>%
  filter(variant_sum > `Total number genotyped` + 0.5) %>%
  dplyr::select(`Record Number`, `Site 1: Country`, `Year study started`, `Year study ended`,
                `Total number genotyped`, variant_sum)
record_variants <- record_variants %>%
  mutate(excess = variant_sum - `Total number genotyped`) %>%
  arrange(desc(excess))
message("CHECK records with more variant counts than Total number genotyped: ", nrow(record_variants),
        " (", round(sum(record_variants$excess)), " extra samples; full list in check_records_overcounted.csv). Largest:")
print(head(record_variants, 10))
write.csv(record_variants, "check_records_overcounted.csv", row.names = FALSE)


# =====================================================================
# 5. Aggregate by country and time period
# =====================================================================
columns_to_sum <- c("A", "B", "C", "D", "F", "G", "H", "J", "K", "L",
                    "N", "O", "P", CRF_cols, "total_CRF", "other_CRF",
                    "URFs", "unspecified_recombinants", "total_recombinants",
                    "Total number genotyped")
columns_to_sum_adjusted <- paste0(columns_to_sum, "_adjusted")

countries_sampled <- HIV_samples %>%
  dplyr::select(
    `Record Number`, no_years, freq_var, year, year_category,
    `Number of records for this publication (ie sites with seperately reported analyses)`,
    starts_with("Site"), `Total number genotyped`, `Study sites description`,
    `Year study started`, `Year study ended`, `Ethnic Group`, `Country of origin`,
    starts_with("Country of origin"), starts_with("If Other, please specify:"),
    all_of(columns_to_sum)
  ) %>%
  mutate(C = as.numeric(C)) %>%
  mutate(across(all_of(columns_to_sum), ~ . / no_years, .names = "{.col}_adjusted")) %>%
  group_by(year_category, `Site 1: Country`) %>%
  dplyr::summarise(across(all_of(columns_to_sum_adjusted), ~ sum(., na.rm = TRUE))) %>%
  arrange(`Site 1: Country`, year_category)

# Denominator: the variant categories that partition the genotyped
# samples. total_CRF, other_CRF and total_recombinants are roll-ups and
# are excluded so they are not double-counted. Their sum,
# total_recorded_genotypes_adj, is the sample size of each country-period
# (FIX 8) and the denominator of proportions.csv and the % columns in
# indices.csv. RMSD and the diversity indices use classified variants
# only (FIX 7, section 6).
subtype_totals <- paste0(
  c("A", "B", "C", "D", "F", "G", "H", "J", "K", "L", "N", "O", "P",
    CRF_cols, "URFs", "unspecified_recombinants"), "_adjusted")

countries_sampled$total_recorded_genotypes_adj <-
  rowSums(countries_sampled[, c(subtype_totals)], na.rm = TRUE)

# CHECK: 174 variant categories (13 subtypes/groups, the CRF columns,
# URFs, unspecified recombinants)
message("CHECK variant categories in the denominator: ", length(subtype_totals))
stopifnot(length(subtype_totals) == 174)

# CHECK: no negative country-period counts, and no samples gained or lost
# between the records and the country-period totals
stopifnot(all(as.matrix(countries_sampled[, columns_to_sum_adjusted]) >= 0))
stopifnot(abs(sum(countries_sampled$total_recorded_genotypes_adj) - records_total_variants) < 1e-6)
message("CHECK total samples (summed variant counts): records = ", records_total_variants,
        "; country-periods = ", round(sum(countries_sampled$total_recorded_genotypes_adj), 1),
        " (country-periods: ", nrow(countries_sampled), ")")

# CHECK (report only): country-periods where the variant counts do not
# add up to Total number genotyped (HIV-2 samples, data entry)
genotyped_mismatch <- countries_sampled %>%
  ungroup() %>%
  mutate(difference = `Total number genotyped_adjusted` - total_recorded_genotypes_adj) %>%
  filter(abs(difference) > 0.5) %>%
  arrange(desc(abs(difference))) %>%
  dplyr::select(`Site 1: Country`, year_category, `Total number genotyped_adjusted`,
                total_recorded_genotypes_adj, difference)
message("CHECK country-periods where variant counts differ from Total number genotyped by > 0.5: ",
        nrow(genotyped_mismatch), " of ", nrow(countries_sampled),
        "; total difference = ", round(sum(genotyped_mismatch$difference), 1), ". Largest:")
print(head(genotyped_mismatch, 10))

write.csv(countries_sampled, "counts_country_subtype_distirbutions.csv")

# Proportions. The roll-ups total_recombinants and total_CRF are converted
# too, for convenience. other_CRF_adjusted and Total number
# genotyped_adjusted stay as counts (the pie charts recompute other_CRF).
# Country-periods with no HIV-1 counts get NaN (0 / 0).
subtype_totals <- paste0(
  c("A", "B", "C", "D", "F", "G", "H", "J", "K", "L", "N", "O", "P",
    CRF_cols, "URFs", "total_recombinants", "total_CRF",
    "unspecified_recombinants"), "_adjusted")

proportions <- countries_sampled %>%
  mutate(across(all_of(subtype_totals), ~ . / total_recorded_genotypes_adj))

write.csv(proportions, "proportions.csv")


# =====================================================================
# 6. Pairwise RMSD
#
# For each pair of countries in each period, the root mean squared
# deviation between their variant proportion vectors. Lower means more
# similar distributions.
#
# FIX 7: RMSD is computed over CLASSIFIED variants only.
# Previously it used 176 columns: the 174 categories plus the roll-ups
# total_recombinants and total_CRF, so a CRF counted three times and the
# maximum RMSD was sqrt(4/176) = 0.1508. Now rmsd_cols holds subtypes
# A-P, the named CRFs and URFs (undefined URFs are inside URFs and count
# as classified): k = 172 columns. Unspecified CRFs and unspecified
# recombinants are left out, and each country's proportions are
# re-scaled to sum to 1 over the classified samples. Country-periods
# with no classified samples get no RMSD. The maximum possible RMSD is
# now sqrt(2/k), checked after the loop.
# proportions.csv above is unchanged, so the pie charts still show
# unspecified recombinants.
# =====================================================================
rmsd_cols <- paste0(
  c("A", "B", "C", "D", "F", "G", "H", "J", "K", "L", "N", "O", "P",
    setdiff(CRF_cols, "Number of unspecified CRFs"), "URFs"), "_adjusted")

classified_props <- countries_sampled %>%
  ungroup() %>%
  mutate(total_classified_adj = rowSums(across(all_of(rmsd_cols)), na.rm = TRUE)) %>%
  filter(total_classified_adj > 0) %>%
  mutate(across(all_of(rmsd_cols), ~ replace_na(.x, 0) / total_classified_adj))

# CHECK: k = 172 classified columns
message("CHECK RMSD columns k = ", length(rmsd_cols), "; country-periods with classified samples: ",
        nrow(classified_props), " of ", nrow(countries_sampled))
stopifnot(length(rmsd_cols) == 172)

calculate_rmsd <- function(df1, df2) {
  if (!all(names(df1) == names(df2))) {
    stop("Dataframes do not have the same structure")
  }
  
  df1_numeric <- df1 %>% ungroup() %>% select(where(is.numeric))
  df2_numeric <- df2 %>% ungroup() %>% select(where(is.numeric))
  
  squared_diff <- (df1_numeric - df2_numeric)^2
  
  sqrt(sum(squared_diff) / ncol(df1_numeric))
}

hiv_data <- classified_props %>%
  select(`Site 1: Country`, year_category, all_of(rmsd_cols)) %>%
  mutate(across(all_of(rmsd_cols), as.numeric))

results <- data.frame(Country_A   = character(),
                      Country_B   = character(),
                      Time_Period = character(),
                      RMSD        = numeric(),
                      stringsAsFactors = FALSE)

time_periods <- unique(hiv_data$year_category)
countries    <- unique(hiv_data$`Site 1: Country`)

for (time_period in time_periods) {
  period_data <- hiv_data %>%
    filter(year_category == time_period) %>%
    select(-year_category)
  
  period_countries <- unique(period_data$`Site 1: Country`)
  
  # a period with only one country has no pairs (without this, 1:(1 - 1)
  # runs as 1:0 and adds a false RMSD of 0 against a missing country)
  if (length(period_countries) < 2) next
  
  for (i in 1:(length(period_countries) - 1)) {
    for (j in (i + 1):length(period_countries)) {
      country_a <- period_countries[i]
      country_b <- period_countries[j]
      
      data_a <- period_data %>%
        filter(`Site 1: Country` == country_a) %>%
        select(-`Site 1: Country`)
      
      data_b <- period_data %>%
        filter(`Site 1: Country` == country_b) %>%
        select(-`Site 1: Country`)
      
      rmsd <- calculate_rmsd(data_a, data_b)
      
      results <- rbind(results, data.frame(Country_A   = country_a,
                                           Country_B   = country_b,
                                           Time_Period = time_period,
                                           RMSD        = rmsd,
                                           stringsAsFactors = FALSE))
    }
  }
}

# check: no RMSD above the maximum possible for k classified columns
message("RMSD over k = ", length(rmsd_cols), " classified variants; max RMSD = ",
        round(max(results$RMSD), 4), "; bound sqrt(2/k) = ", round(sqrt(2 / length(rmsd_cols)), 4))
stopifnot(max(results$RMSD) <= sqrt(2 / length(rmsd_cols)) + 1e-12)
# CHECK: no missing or negative RMSD, and one row for every country pair
expected_pairs <- hiv_data %>% count(year_category) %>% summarise(pairs = sum(choose(n, 2)))
message("CHECK RMSD rows: ", nrow(results), " (expected ", expected_pairs$pairs, ")")
stopifnot(all(!is.na(results$RMSD)), all(results$RMSD >= 0),
          nrow(results) == expected_pairs$pairs,
          !any(duplicated(results[, c("Country_A", "Country_B", "Time_Period")])))

write.csv(results, "RMSDcountries.csv")


# =====================================================================
# 7. Country-level diversity indices
# =====================================================================
# FIX 7 (continued): Shannon, Simpson and the variant / CRF counts use
# the same classified columns as RMSD (rmsd_cols), so unspecified CRFs
# and unspecified recombinants no longer count as a "variant" or a
# "CRF". The % unclassified is reported separately below.
counts_data <- countries_sampled %>%
  ungroup() %>%
  dplyr::select(all_of(rmsd_cols)) %>%
  as.matrix()

simpson_indices <- diversity(counts_data, index = "simpson")
shannon_indices <- diversity(counts_data, index = "shannon")

indices <- countries_sampled %>%
  ungroup() %>%
  cbind(simpson_indices, shannon_indices)

# vegan's "simpson" is already 1 - sum(p^2) (Gini-Simpson: the chance
# two random samples are different variants). That is what Table 1
# reports. Note: inv_simpson = 1 - that is
# sum(p^2), Simpson's original D (name kept as inv_simpson).
indices <- indices %>% mutate(inv_simpson = 1 - simpson_indices)

indices$count_subtypes <- rowSums(indices[rmsd_cols] > 0)

CRF_cols_adj <- paste0(setdiff(CRF_cols, "Number of unspecified CRFs"), "_adjusted")
indices$count_CRFs <- rowSums(indices[CRF_cols_adj] > 0)

indices <- indices %>%
  # FIX 8: divide by the summed variant counts, as proportions.csv does,
  # not by the raw Total number genotyped field
  mutate(proportion_recombinants = (total_recombinants_adjusted / total_recorded_genotypes_adj) * 100,
         proportion_CRFs         = (total_CRF_adjusted          / total_recorded_genotypes_adj) * 100,
         proportion_URFs         = (URFs_adjusted               / total_recorded_genotypes_adj) * 100,
         # % of recorded samples that are unclassified (unspecified CRFs
         # + unspecified recombinants), for the tables
         proportion_unclassified = ((`Number of unspecified CRFs_adjusted` + unspecified_recombinants_adjusted) /
                                      total_recorded_genotypes_adj) * 100,
         total_classified_adj    = rowSums(across(all_of(rmsd_cols)), na.rm = TRUE))

# FIX 6: drop country-periods with no HIV-1 variant counts. Cabo Verde
# 2015-2019 holds only an HIV-2 record (record 595), so it had 11.5
# "samples", 0 variants, Shannon 0 and Simpson 1 (vegan's value for an
# all-zero row). It never had an RMSD, but it entered the all-country
# diversity table.
# FIX 7: likewise drop country-periods with no CLASSIFIED samples
indices <- indices %>%
  filter(total_recorded_genotypes_adj > 0, total_classified_adj > 0)

# CHECK: diversity indices in range. Shannon is between 0 and
# log(number of variants); Simpson (1 - sum p^2) is between 0 and 1
message("CHECK indices: Shannon ", round(min(indices$shannon_indices), 3), " to ",
        round(max(indices$shannon_indices), 3), "; Simpson ", round(min(indices$simpson_indices), 3),
        " to ", round(max(indices$simpson_indices), 3), "; % unclassified ",
        round(min(indices$proportion_unclassified), 2), " to ", round(max(indices$proportion_unclassified), 2))
stopifnot(all(!is.na(indices$shannon_indices)), all(!is.na(indices$simpson_indices)),
          all(indices$shannon_indices >= 0),
          all(indices$shannon_indices <= log(indices$count_subtypes) + 1e-9),
          all(indices$simpson_indices >= 0), all(indices$simpson_indices < 1),
          all(indices$proportion_unclassified >= 0), all(indices$proportion_unclassified <= 100))

indices <- indices %>%
  dplyr::select(year_category, `Site 1: Country`, simpson_indices, shannon_indices,
                inv_simpson, count_subtypes, count_CRFs, proportion_recombinants,
                proportion_CRFs, proportion_URFs, proportion_unclassified,
                total_recorded_genotypes_adj)   # FIX 8: sample size = summed variant counts

write.csv(indices, "indices.csv")