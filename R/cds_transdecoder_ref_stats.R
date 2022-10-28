# how many transcripts.
# How many are novel?
# How many have CDS
# How many are NMDj
# How many with full lenght evidence?
# and sets
# Waterfall plot
library(dplyr)
library(GenomicFeatures)
library(ggplot2)
library(ROCR)
library(caret)

compared <- rtracklayer::import("phase2/compared/cds.combined.gtf")
compared <- subset(compared, type == 'transcript')
trd <- rtracklayer::import(
  "phase2/transdecoder/longest_orfs.cds.best_candidates.genome.gff3")

# Total transdecoder CDS
length(unique(trd[trd$type == 'CDS']$ID))

trd <- subset(trd, type == "mRNA")
trd <- trd[grepl(x=trd$Name, pattern = 'ORF type:complete') ,]
trd$score <- gsub(trd$Name, pattern = ".*,score=([0-9.]+)$", replacement = "\\1")

sum(trd$type == 'mRNA')
#n_distinct(trd[trd$type == 'CDS']$ID)
length(compared)
sort(table(compared$class_code), decreasing = T)
# 245,447 transcripts with CDS
# 170,660 of which are complete CDS
# 193,108 comparable with ref
# =     j     u     c     i      k     x     e     p     o     n      m     y     s 
# 64447 53176 32873 18443 10129  5751  2638  1614  1341  1180  1063   299   136   18 

compared_df <- compared %>% 
  mcols() %>% 
  as_tibble() %>% 
  select(-c(score, phase, exon_number, tss_id))

trd_df <- trd %>% 
  mcols() %>% 
  as_tibble() %>% 
  dplyr::select(ID, score) %>% 
  dplyr::rename(oId = ID)

x <- left_join(trd_df, compared_df) %>% 
  mutate(score = as.numeric(score))

x %>% 
  filter(class_code %in% c('=', 'j', 'u', 'c', NA)) %>% 
  ggplot(aes(x=score)) + 
  lims(x=c(-5, 1000)) +
  geom_histogram(binwidth = 4) +
  facet_wrap(~class_code, nrow = 5) 

cuttoff <- x %>% 
  filter(class_code == '=') %>% 
  pull(score) %>% 
  quantile(., 0.05, na.rm=TRUE)
# cutoff for attaining 95% of the identical matches is 23.45 
trd_score_density <- density(x$score, na.rm=TRUE)
trd_score_density <- as.data.frame(trd_score_density[c("x", "y")])
trd_score_density <- subset(trd_score_density, x > cuttoff)

x %>% 
  filter(class_code %in% c('=', 'j', 'u', 'c', NA)) %>% 
  ggplot(aes(x=score, fill=class_code, color=class_code)) + 
  geom_density(alpha=.30) +
  lims(x=c(-10, 100)) + 
  geom_vline(
    xintercept = cuttoff, 
    linetype="dotted", 
    color = "red", 
    size=1) +
  geom_text(
    aes(
      x=cuttoff, 
      label=stringr::str_glue("5% quantile ({cuttoff})"), 
      y=0.06), 
    colour="red", 
    check_overlap = TRUE)
# interestingly class_code of U and NA have a should < cutoff
# other are stable > cutoff
  
x2 <- x %>% filter(class_code %in% c("=", "u"))
labels = ifelse(x2$class_code == "=", 1, -1)
#labels[is.na(labels)] <- -1
predictions = scale(x2$score)
predictions[is.na(predictions)] <- min(predictions, na.rm = T)
pred <- prediction(predictions, labels)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf, colorize=TRUE)
abline(a = 0, b = 1)

auc_ROCR <- performance(pred, measure = "auc")
auc_ROCR <- auc_ROCR@y.values[[1]]

data.frame(x=perf@x.values[[1]],y=perf@y.values[[1]]) %>% 
  ggplot(aes(x=x,y=y))  +
   geom_path(size=1) +
   geom_segment(aes(x=0,y=0,xend=1,yend=1),colour="black",linetype= 2) +
   annotate("text", x=0.90, y= 0, label=paste(sep = "", "AUC = ", round(auc_ROCR,2)),colour="black") +
   scale_x_continuous(name= "False positive rate") +
   scale_y_continuous(name= "True positive rate") + 
  theme_bw()

xtab <- table(
  factor(ifelse(x2$score > cuttoff, 1, -1), levels = c(1, -1)), 
  factor(labels, levels = c(1, -1)))
confusionMatrix(xtab)

ggplot(x2, aes(x = score)) + 
  geom_histogram(binwidth = .05) + 
  # facet_wrap(~obs) + 
  xlab("Probability of Class #1")

# maybe check this: https://cran.r-project.org/web/packages/OptimalCutpoints/

# second test
x3 <- x 
labels = ifelse(x3$class_code == "=", 1, -1)
labels[is.na(labels)] <- -1
predictions = scale(x3$score)
predictions[is.na(predictions)] <- min(predictions, na.rm = T)
pred <- prediction(predictions, labels)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf, colorize=TRUE)
abline(a = 0, b = 1)

auc_ROCR <- performance(pred, measure = "auc")
auc_ROCR <- auc_ROCR@y.values[[1]]

data.frame(x=perf@x.values[[1]],y=perf@y.values[[1]]) %>% 
  ggplot(aes(x=x,y=y))  +
  geom_path(size=1) +
  geom_segment(aes(x=0,y=0,xend=1,yend=1),colour="black",linetype= 2) +
  annotate("text", x=0.90, y= 0, label=paste(sep = "", "AUC = ", round(auc_ROCR,2)),colour="black") +
  scale_x_continuous(name= "False positive rate") +
  scale_y_continuous(name= "True positive rate") + 
  theme_bw()

xtab <- table(
  factor(ifelse(x3$score > cuttoff, 1, -1), levels = c(1, -1)), 
  factor(labels, levels = c(1, -1)))
confusionMatrix(xtab)
