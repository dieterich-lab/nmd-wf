
metadata <- read_csv(
  "phase2/config/metadata_w_files.csv", 
  col_types = cols(
    ...1 = col_double(),
    CCG_Sample_ID = col_character(),
    Condition = col_character(),
    VB_Internal_Sample_ID = col_character(),
    Cell_line = col_character(),
    Knockdown = col_character(),
    Replicate = col_double(),
    Relevant = col_logical(),
    `Spike-Ins` = col_character(),
    Library_prep = col_character(),
    Sequencing = col_character(),
    platform = col_character(),
    bams = col_character(),
    raw_reads = col_character()))
  
metadata %>%
  select(-c(...1)) %>%
  mutate(
    Knockdown = str_replace(Knockdown, '( +)', ''),
    cellline = word(Cell_line, 1),
    Condition = str_replace(Condition, '( +)', ''),
    Knockout = str_extract(Condition, "SMG7"),
    Knockout = ifelse(Knockout == 'SMG7', "_SMG7-KO", ""),
    clone = str_extract(Cell_line,  "(?<=\\().+?(?=\\))"),
    group = str_glue("{cellline}{Knockout}_{Knockdown}-KD_{clone}", .na=''),
    group = str_replace(group, '(_$)', ''))  %>%
  select(CCG_Sample_ID, VB_Internal_Sample_ID, group, everything()) %>% 
  write_csv("phase2/config/metadata_w_files.csv")
