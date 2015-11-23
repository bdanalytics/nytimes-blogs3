# NYTimes:Blogs: Popular classification:: NYTBlogs3_feat_PubDate
bdanalytics  

**  **    
**Date: (Mon) Nov 23, 2015**    

# Introduction:  

Data: 
Source: 
    Training:   https://www.kaggle.com/c/15-071x-the-analytics-edge-competition-spring-2015/download/NYTimesBlogTrain.csv  
    New:        https://www.kaggle.com/c/15-071x-the-analytics-edge-competition-spring-2015/download/NYTimesBlogTest.csv  
Time period: 



# Synopsis:

Based on analysis utilizing <> techniques, <conclusion heading>:  

Summary of key steps & error improvement stats:

### Prediction Accuracy Enhancement Options:
- transform.data chunk:
    - derive features from multiple features
    
- manage.missing.data chunk:
    - Not fill missing vars
    - Fill missing numerics with a different algorithm
    - Fill missing chars with data based on clusters 
    
- extract.features chunk:
    - Text variables: move to date extraction chunk ???
        - Mine acronyms
        - Mine places

- Review set_global_options chunk after features are finalized

### ![](<filename>.png)

## Potential next steps include:
- Organization:
    - Categorize by chunk
    - Priority criteria:
        0. Ease of change
        1. Impacts report
        2. Cleans innards
        3. Bug report
        
- all chunks:
    - at chunk-end rm(!glb_<var>)
    
- manage.missing.data chunk:
    - cleaner way to manage re-splitting of training vs. new entity

- extract.features chunk:
    - Add n-grams for glbFeatsText
        - "RTextTools", "tau", "RWeka", and "textcat" packages
    - Convert user-specified mutate code to config specs
    
- fit.models chunk:
    - Prediction accuracy scatter graph:
    -   Add tiles (raw vs. PCA)
    -   Use shiny for drop-down of "important" features
    -   Use plot.ly for interactive plots ?
    
    - Change .fit suffix of model metrics to .mdl if it's data independent (e.g. AIC, Adj.R.Squared - is it truly data independent ?, etc.)
    - create a custom model for rpart that has minbucket as a tuning parameter
    - varImp for randomForest crashes in caret version:6.0.41 -> submit bug report

- Probability handling for multinomials vs. desired binomial outcome
-   ROCR currently supports only evaluation of binary classification tasks (version 1.0.7)
-   extensions toward multiclass classification are scheduled for the next release

- Skip trControl.method="cv" for dummy classifier ?
- Add custom model to caret for a dummy (baseline) classifier (binomial & multinomial) that generates proba/outcomes which mimics the freq distribution of glb_rsp_var values; Right now glb_dmy_glm_mdl always generates most frequent outcome in training data
- glm_dmy_mdl should use the same method as glm_sel_mdl until custom dummy classifer is implemented

- fit.all.training chunk:
    - myplot_prediction_classification: displays 'x' instead of '+' when there are no prediction errors 
- Compare glb_sel_mdl vs. glb_fin_mdl:
    - varImp
    - Prediction differences (shd be minimal ?)

- Move glb_analytics_diag_plots to mydsutils.R: (+) Easier to debug (-) Too many glb vars used
- Add print(ggplot.petrinet(glb_analytics_pn) + coord_flip()) at the end of every major chunk
- Parameterize glb_analytics_pn
- Move glb_impute_missing_data to mydsutils.R: (-) Too many glb vars used; glb_<>_df reassigned
- Do non-glm methods handle interaction terms ?
- f-score computation for classifiers should be summation across outcomes (not just the desired one ?)
- Add accuracy computation to glb_dmy_mdl in predict.data.new chunk
- Why does splitting fit.data.training.all chunk into separate chunks add an overhead of ~30 secs ? It's not rbind b/c other chunks have lower elapsed time. Is it the number of plots ?
- Incorporate code chunks in print_sessionInfo
- Test against 
    - projects in github.com/bdanalytics
    - lectures in jhu-datascience track

# Analysis: 

```r
rm(list = ls())
set.seed(12345)
options(stringsAsFactors = FALSE)
source("~/Dropbox/datascience/R/myscript.R")
source("~/Dropbox/datascience/R/mydsutils.R")
```

```
## Loading required package: caret
## Loading required package: lattice
## Loading required package: ggplot2
```

```r
source("~/Dropbox/datascience/R/myplot.R")
source("~/Dropbox/datascience/R/mypetrinet.R")
source("~/Dropbox/datascience/R/myplclust.R")
source("~/Dropbox/datascience/R/mytm.R")
# Gather all package requirements here
suppressPackageStartupMessages(require(doMC))
registerDoMC(6) # # of cores on machine - 2
suppressPackageStartupMessages(require(caret))
#source("dbgcaret.R")
#packageVersion("snow")
#require(sos); findFn("cosine", maxPages=2, sortby="MaxScore")

# Analysis control global variables
# Inputs
glb_trnng_url <- "https://www.kaggle.com/c/15-071x-the-analytics-edge-competition-spring-2015/download/NYTimesBlogTrain.csv"
glb_newdt_url <- "https://www.kaggle.com/c/15-071x-the-analytics-edge-competition-spring-2015/download/NYTimesBlogTest.csv"
glbInpMerge <- NULL #: default
#     list(fnames = c("<fname1>", "<fname2>")) # files will be concatenated

glb_is_separate_newobs_dataset <- TRUE    # or TRUE
    glb_split_entity_newobs_datasets <- FALSE   # select from c(FALSE, TRUE)
    glb_split_newdata_method <- NULL # select from c(NULL, "condition", "sample", "copy")
    glb_split_newdata_condition <- NULL # or "is.na(<var>)"; "<var> <condition_operator> <value>"
    glb_split_newdata_size_ratio <- 0.3               # > 0 & < 1
    glb_split_sample.seed <- 123               # or any integer

glbObsDropCondition <- NULL # : default
#            "<condition>" # use | & ; NOT || &&
#parse(text=glbObsDropCondition)
#subset(glbObsAll, .grpid %in% c(31))
    
glb_obs_repartition_train_condition <- NULL # : default
#    "<condition>" 

glb_max_fitobs <- NULL # or any integer
                         
glb_is_regression <- FALSE; glb_is_classification <- !glb_is_regression; 
    glb_is_binomial <- TRUE # or TRUE or FALSE

glb_rsp_var_raw <- "Popular"

# for classification, the response variable has to be a factor
glb_rsp_var <- "Pplr.fctr" # glb_rsp_var_raw # or "Pplr.fctr"

# if the response factor is based on numbers/logicals e.g (0/1 OR TRUE/FALSE vs. "A"/"B"), 
#   or contains spaces (e.g. "Not in Labor Force")
#   caret predict(..., type="prob") crashes
glb_map_rsp_raw_to_var <- #NULL 
function(raw) {
#     return(raw ^ 0.5)
#     return(log(1 + raw))
#     return(log10(raw)) 
#     return(exp(-raw / 2))
    ret_vals <- rep_len(NA, length(raw)); ret_vals[!is.na(raw)] <- ifelse(raw[!is.na(raw)] == 1, "Y", "N"); return(relevel(as.factor(ret_vals), ref="N"))
#     #as.factor(paste0("B", raw))
#     #as.factor(gsub(" ", "\\.", raw))    
    }
# glb_map_rsp_raw_to_var(tst <- c(NA, 0, 1)) 
# glb_map_rsp_raw_to_var(tst <- c(NA, 0, 2.99, 280.50, 1000.00)) 

glb_map_rsp_var_to_raw <- #NULL 
function(var) {
#     return(var ^ 2.0)
#     return(exp(var))
#     return(10 ^ var) 
#     return(-log(var) * 2)
    as.numeric(var) - 1
#     gsub("\\.", " ", levels(var)[as.numeric(var)])
#     c("<=50K", " >50K")[as.numeric(var)]
#     c(FALSE, TRUE)[as.numeric(var)]
}
# glb_map_rsp_var_to_raw(glb_map_rsp_raw_to_var(tst))

if ((glb_rsp_var != glb_rsp_var_raw) && is.null(glb_map_rsp_raw_to_var))
    stop("glb_map_rsp_raw_to_var function expected")

# List info gathered for various columns
# <col_name>:   <description>; <notes>
# NewsDesk = the New York Times desk that produced the story (Business, Culture, Foreign, etc.)
# SectionName = the section the article appeared in (Opinion, Arts, Technology, etc.)
# SubsectionName = the subsection the article appeared in (Education, Small Business, Room for Debate, etc.)
# Headline = the title of the article
# Snippet = a small portion of the article text
# Abstract = a summary of the blog article, written by the New York Times
# WordCount = the number of words in the article
# PubDate = the publication date, in the format "Year-Month-Day Hour:Minute:Second"
# UniqueID = a unique identifier for each article

# If multiple vars are parts of id, consider concatenating them to create one id var
# If glb_id_var == NULL, ".rownames <- row.names()" is the default

# User-specified exclusions
glbFeatsExclude <- c(NULL
#   Feats that shd be excluded due to known causation by prediction variable
# , "<feat1", "<feat2>"
#   Feats that are linear combinations (alias in glm)
#   Feature-engineering phase -> start by excluding all features except id & category & work each one in
    , "NewsDesk", "SectionName", "SubsectionName"
    , "WordCount", "PubDate" 
    # Feature Engineering done with prior features
    , "Headline", "Snippet", "Abstract"
                    ) 
if (glb_rsp_var_raw != glb_rsp_var)
    glbFeatsExclude <- union(glbFeatsExclude, glb_rsp_var_raw)                    

glbFeatsInteractionOnly <- list()
#glbFeatsInteractionOnly[["<child_feat>"]] <- "<parent_feat>"

# currently does not handle more than 1 column; consider concatenating multiple columns
glb_id_var <- "UniqueID" # choose from c(NULL : default, "<id_feat>") 
glbFeatsCategory <- "NDSS.my.fctr" # choose from c(NULL : default, "<category>")

glb_drop_vars <- c(NULL
                # , "<feat1>", "<feat2>"
                )

glb_map_vars <- NULL # or c("<var1>", "<var2>")
glb_map_urls <- list();
# glb_map_urls[["<var1>"]] <- "<var1.url>"

glb_assign_pairs_lst <- NULL; 
# glb_assign_pairs_lst[["<var1>"]] <- list(from=c(NA),
#                                            to=c("NA.my"))
glb_assign_vars <- names(glb_assign_pairs_lst)

# Derived features; Use this mechanism to cleanse data ??? Cons: Data duplication ???
glbFeatsDerive <- list();

# glbFeatsDerive[["<feat.my.sfx>"]] <- list(
#     mapfn = function(<arg1>, <arg2>) { return(function(<arg1>, <arg2>)) } 
#   , args = c("<arg1>", "<arg2>"))

    # character
#     mapfn = function(Week) { return(substr(Week, 1, 10)) }

#     mapfn = function(descriptor) { return(plyr::revalue(descriptor, c(
#         "ABANDONED BUILDING"  = "OTHER",
#         "**"                  = "**"
#                                           ))) }
glbFeatsDerive[["NDSS.my"]] <- list(
    mapfn = function(NewsDesk, SectionName, SubsectionName) { 
        descriptor <- 
            gsub(" ", "", paste(NewsDesk, SectionName, SubsectionName, sep = "#"))
        return(gsub("[a|e|i|o|u]", "", plyr::revalue(descriptor, c(NULL
            , "#BusinessDay#Dealbook"       = "Business#BusinessDay#Dealbook"
            , "#BusinessDay#SmallBusiness"  = "Business#BusinessDay#SmallBusiness"
            , "#Crosswords/Games#"          = "Business#Crosswords/Games#"
            , "#Open#"                      = "Business#Technology#"            
            , "#Technology#"                = "Business#Technology#"
            , "Business##"                  = "Business#Technology#"            
            
            , "#Arts#"                      = "Culture#Arts#"            
            
            , "Foreign##"           = "Foreign#World#"
            , "#World#AsiaPacific"  = "Foreign#World#AsiaPacific"  
            
            , "#N.Y./Region#"       = "Metro#N.Y./Region#"
            
            , "#Opinion#"           = "OpEd#Opinion#"
            , "OpEd##"              = "OpEd#Opinion#"
            
            , "#Health#"            = "Science#Health#"
            , "Science##"           = "Science#Health#"            
            , "Styles#Health#"      = "Science#Health#"
            
            , "Styles##"                   = "Styles##Fashion"
            , "Styles#Style#Fashion&Style" = "Styles##Fashion"
            
            , "#Travel#"                = "Travel#Travel#"            
            
            , "Magazine#Magazine#"      = "myOther"            
            , "National##"              = "myOther"
            , "National#U.S.#Politics"  = "myOther"        
            , "Sports##"                = "myOther"
            , "Sports#Sports#"          = "myOther"
            , "#U.S.#"                  = "myOther"         
                                          ))))
        }
    , args = c("NewsDesk", "SectionName", "SubsectionName"))    

#     mapfn = function(description) { mod_raw <- description;
    # This is here because it does not work if it's in txt_map_filename
#         mod_raw <- gsub(paste0(c("\n", "\211", "\235", "\317", "\333"), collapse = "|"), " ", mod_raw)
    # Don't parse for "." because of ".com"; use customized gsub for that text
#         mod_raw <- gsub("(\\w)(!|\\*|,|-|/)(\\w)", "\\1\\2 \\3", mod_raw);
#         return(mod_raw) }
#print(mod_raw <- grep("&#034;", glbObsAll[, txt_var], value = TRUE)) 
#print(mod_raw <- glbObsAll[c(88,187,280,1040,1098), txt_var])
#print(mod_raw <- glbObsAll[sel_obs(list(descr.my.contains="\\bdoes( +)not\\b")), glbFeatsText])
#print(mod_raw <- glbObsAll[sel_obs(list(descr.my.contains="\\bipad [[:digit:]]\\b")), glbFeatsText][01:10])
#print(mod_raw <- glbObsAll[sel_obs(list(descr.my.contains="pad mini")), glbFeatsText][11:20])
#print(mod_raw <- glbObsAll[sel_obs(list(descr.my.contains="pad mini")), glbFeatsText][21:30])
#print(mod_raw <- glbObsAll[sel_obs(list(descr.my.contains="pad mini")), glbFeatsText][31:40])
#glbObsAll[which(glb_post_stop_words_terms_mtrx_lst[[txt_var]][, subset(glb_post_stop_words_terms_df_lst[[txt_var]], term %in% c("conditionminimal"))$pos] > 0), "description"]

    # numeric
# Create feature based on record position/id in data   
# glbFeatsDerive[["dummy.my"]] <- list(
#     mapfn = function(UniqueID) { return(UniqueID) }       
#     , args = c("UniqueID"))    
    
# Add logs of numerics that are not distributed normally
#   Derive & keep multiple transformations of the same feature, if normality is hard to achieve with just one transformation
#   Right skew: logp1; sqrt; ^ 1/3; logp1(logp1); log10; exp(-<feat>/constant)
glbFeatsDerive[["WordCount.log1p"]] <- list(
    mapfn = function(WordCount) { return(log1p(WordCount)) } 
  , args = c("WordCount"))
glbFeatsDerive[["WordCount.root2"]] <- list(
    mapfn = function(WordCount) { return(WordCount ^ (1/2)) } 
  , args = c("WordCount"))
glbFeatsDerive[["WordCount.nexp"]] <- list(
    mapfn = function(WordCount) { return(exp(-WordCount)) } 
  , args = c("WordCount"))
#print(summary(glbObsAll$WordCount))
#print(summary(mapfn(glbObsAll$WordCount)))
    
#     mapfn = function(Rasmussen) { return(ifelse(sign(Rasmussen) >= 0, 1, 0)) } 
#     mapfn = function(startprice) { return(startprice ^ (1/2)) }       
#     mapfn = function(startprice) { return(log(startprice)) }   
#     mapfn = function(startprice) { return(exp(-startprice / 20)) }
#     mapfn = function(startprice) { return(scale(log(startprice))) }     
#     mapfn = function(startprice) { return(sign(sprice.predict.diff) * (abs(sprice.predict.diff) ^ (1/10))) }        

    # factor      
#     mapfn = function(PropR) { return(as.factor(ifelse(PropR >= 0.5, "Y", "N"))) }
#     mapfn = function(productline, description) { as.factor(gsub(" ", "", productline)) }
#     mapfn = function(purpose) { return(relevel(as.factor(purpose), ref="all_other")) }
#     mapfn = function(raw) { tfr_raw <- as.character(cut(raw, 5)); 
#                             tfr_raw[is.na(tfr_raw)] <- "NA.my";
#                             return(as.factor(tfr_raw)) }
#     mapfn = function(startprice.log10) { return(cut(startprice.log10, 3)) }
#     mapfn = function(startprice.log10) { return(cut(sprice.predict.diff, c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000))) }    

#     , args = c("<arg1>"))
    
    # multiple args    
#     mapfn = function(PTS, oppPTS) { return(PTS - oppPTS) }
#     mapfn = function(startprice.log10.predict, startprice) {
#                  return(spdiff <- (10 ^ startprice.log10.predict) - startprice) } 
#     mapfn = function(productline, description) { as.factor(
#         paste(gsub(" ", "", productline), as.numeric(nchar(description) > 0), sep = "*")) }

# # If glbObsAll is not sorted in the desired manner
#     mapfn=function(Week) { return(coredata(lag(zoo(orderBy(~Week, glbObsAll)$ILI), -2, na.pad=TRUE))) }
#     mapfn=function(ILI) { return(coredata(lag(zoo(ILI), -2, na.pad=TRUE))) }
#     mapfn=function(ILI.2.lag) { return(log(ILI.2.lag)) }

# glbFeatsDerive[["<var1>"]] <- glbFeatsDerive[["<var2>"]]

glb_derive_vars <- names(glbFeatsDerive)

# tst <- "descr.my"; args_lst <- NULL; for (arg in glbFeatsDerive[[tst]]$args) args_lst[[arg]] <- glbObsAll[, arg]; print(head(args_lst[[arg]])); print(head(drv_vals <- do.call(glbFeatsDerive[[tst]]$mapfn, args_lst))); 
# print(which_ix <- which(args_lst[[arg]] == 0.75)); print(drv_vals[which_ix]); 

glbFeatsDateTime <- list()
glbFeatsDateTime[["PubDate"]] <- 
    c(format = "%Y-%m-%d %H:%M:%S", timezone = "America/New_York", impute.na = FALSE, 
      last.ctg = FALSE, poly.ctg = FALSE)

glbFeatsPrice <- NULL # or c("<price_var>")

glbFeatsText <- NULL # c("<txt_var>")   # NULL # 
Sys.setlocale("LC_ALL", "C") # For english
```

```
## [1] "C/C/C/C/C/en_US.UTF-8"
```

```r
# Text Processing Step: custom modifications not present in txt_munge -> use glbFeatsDerive
# Text Processing Step: universal modifications
glb_txt_munge_filenames_pfx <- "NYTBlogs3_mytxt_"

# Text Processing Step: tolower
# Text Processing Step: myreplacePunctuation
# Text Processing Step: removeWords
glb_txt_stop_words <- list()
# Remember to use unstemmed words
if (!is.null(glbFeatsText)) {
    require(tm)

    glb_txt_stop_words[["<txt_var>"]] <- sort(c(NULL    

        # Remove any words from stopwords            
#         , setdiff(myreplacePunctuation(stopwords("english")), c("<keep_wrd1>", <keep_wrd2>"))
                                
        # cor.y.train == NA
#         ,unlist(strsplit(paste(c(NULL
#           ,"<comma-separated-terms>"
#         ), collapse=",")

        # freq == 1; keep c("<comma-separated-terms-to-keep>")
            # ,<comma-separated-terms>

        # chisq.pval high (e.g. == 1); keep c("<comma-separated-terms-to-keep>")
            # ,<comma-separated-terms>

        # nzv.freqRatio high (e.g. >= glb_nzv_freqCut); keep c("<comma-separated-terms-to-keep>")
            # ,<comma-separated-terms>        
                                            ))
}
#orderBy(~term, glb_post_stem_words_terms_df_lst[[txt_var]][grep("^2", glb_post_stem_words_terms_df_lst[[txt_var]]$term), ])
#glbObsAll[glb_post_stem_words_terms_mtrx_lst[[txt_var]][, 6] > 0, glbFeatsText]

# To identify terms with a specific freq
#paste0(sort(subset(glb_post_stop_words_terms_df_lst[[txt_var]], freq == 1)$term), collapse = ",")
#paste0(sort(subset(glb_post_stem_words_terms_df_lst[[txt_var]], freq <= 2)$term), collapse = ",")

# To identify terms with a specific freq & 
#   are not stemmed together later OR is value of color.fctr (e.g. gold)
#paste0(sort(subset(glb_post_stop_words_terms_df_lst[[txt_var]], (freq == 1) & !(term %in% c("blacked","blemish","blocked","blocks","buying","cables","careful","carefully","changed","changing","chargers","cleanly","cleared","connect","connects","connected","contains","cosmetics","default","defaulting","defective","definitely","describe","described","devices","displays","drop","drops","engravement","excellant","excellently","feels","fix","flawlessly","frame","framing","gentle","gold","guarantee","guarantees","handled","handling","having","install","iphone","iphones","keeped","keeps","known","lights","line","lining","liquid","liquidation","looking","lots","manuals","manufacture","minis","most","mostly","network","networks","noted","opening","operated","performance","performs","person","personalized","photograph","physically","placed","places","powering","pre","previously","products","protection","purchasing","returned","rotate","rotation","running","sales","second","seconds","shipped","shuts","sides","skin","skinned","sticker","storing","thats","theres","touching","unusable","update","updates","upgrade","weeks","wrapped","verified","verify") ))$term), collapse = ",")

#print(subset(glb_post_stem_words_terms_df_lst[[txt_var]], (freq <= 2)))
#glbObsAll[which(terms_mtrx[, 229] > 0), glbFeatsText]

# To identify terms with cor.y == NA
#orderBy(~-freq+term, subset(glb_post_stop_words_terms_df_lst[[txt_var]], is.na(cor.y)))
#paste(sort(subset(glb_post_stop_words_terms_df_lst[[txt_var]], is.na(cor.y))[, "term"]), collapse=",")
#orderBy(~-freq+term, subset(glb_post_stem_words_terms_df_lst[[txt_var]], is.na(cor.y)))

# To identify terms with low cor.y.abs
#head(orderBy(~cor.y.abs+freq+term, subset(glb_post_stem_words_terms_df_lst[[txt_var]], !is.na(cor.y))), 5)

# To identify terms with high chisq.pval
#subset(glb_post_stem_words_terms_df_lst[[txt_var]], chisq.pval > 0.99)
#paste0(sort(subset(glb_post_stem_words_terms_df_lst[[txt_var]], (chisq.pval > 0.99) & (freq <= 10))$term), collapse=",")
#paste0(sort(subset(glb_post_stem_words_terms_df_lst[[txt_var]], (chisq.pval > 0.9))$term), collapse=",")
#head(orderBy(~-chisq.pval+freq+term, glb_post_stem_words_terms_df_lst[[txt_var]]), 5)
#glbObsAll[glb_post_stem_words_terms_mtrx_lst[[txt_var]][, 68] > 0, glbFeatsText]
#orderBy(~term, glb_post_stem_words_terms_df_lst[[txt_var]][grep("^m", glb_post_stem_words_terms_df_lst[[txt_var]]$term), ])

# To identify terms with high nzv.freqRatio
#summary(glb_post_stem_words_terms_df_lst[[txt_var]]$nzv.freqRatio)
#paste0(sort(setdiff(subset(glb_post_stem_words_terms_df_lst[[txt_var]], (nzv.freqRatio >= glb_nzv_freqCut) & (freq < 10) & (chisq.pval >= 0.05))$term, c( "128gb","3g","4g","gold","ipad1","ipad3","ipad4","ipadair2","ipadmini2","manufactur","spacegray","sprint","tmobil","verizon","wifion"))), collapse=",")

# To identify obs with a txt term
#tail(orderBy(~-freq+term, glb_post_stop_words_terms_df_lst[[txt_var]]), 20)
#mydspObs(list(descr.my.contains="non"), cols=c("color", "carrier", "cellular", "storage"))
#grep("ever", dimnames(terms_stop_mtrx)$Terms)
#which(terms_stop_mtrx[, grep("ipad", dimnames(terms_stop_mtrx)$Terms)] > 0)
#glbObsAll[which(terms_stop_mtrx[, grep("16", dimnames(terms_stop_mtrx)$Terms)[1]] > 0), c(glbFeatsCategory, "storage", txt_var)]

# To identify whether terms shd be synonyms
#orderBy(~term, glb_post_stop_words_terms_df_lst[[txt_var]][grep("^moder", glb_post_stop_words_terms_df_lst[[txt_var]]$term), ])
# term_row_df <- glb_post_stop_words_terms_df_lst[[txt_var]][grep("^came$", glb_post_stop_words_terms_df_lst[[txt_var]]$term), ]
# 
# cor(glb_post_stop_words_terms_mtrx_lst[[txt_var]][glbObsAll$.lcn == "Fit", term_row_df$pos], glbObsTrn[, glb_rsp_var], use="pairwise.complete.obs")

# To identify which stopped words are "close" to a txt term
#sort(cluster_vars)

# Text Processing Step: stemDocument
# To identify stemmed txt terms
#glb_post_stop_words_terms_df_lst[[txt_var]][grep("condit", glb_post_stop_words_terms_df_lst[[txt_var]]$term), ]
#orderBy(~term, glb_post_stem_words_terms_df_lst[[txt_var]][grep("^con", glb_post_stem_words_terms_df_lst[[txt_var]]$term), ])
#glbObsAll[which(terms_stem_mtrx[, grep("use", dimnames(terms_stem_mtrx)$Terms)[[1]]] > 0), c(glb_id_var, "productline", txt_var)]
#glbObsAll[which(TfIdf_stem_mtrx[, 191] > 0), c(glb_id_var, glbFeatsCategory, txt_var)]
#which(glbObsAll$UniqueID %in% c(11915, 11926, 12198))

# Text Processing Step: mycombineSynonyms
#   To identify which terms are associated with not -> combine "could not" & "couldn't"
#findAssocs(glb_full_DTM_lst[[txt_var]], "not", 0.05)
#   To identify which synonyms should be combined
#orderBy(~term, glb_post_stem_words_terms_df_lst[[txt_var]][grep("^c", glb_post_stem_words_terms_df_lst[[txt_var]]$term), ])
chk_comb_cor <- function(syn_lst) {
#     cor(terms_stem_mtrx[glbObsAll$.src == "Train", grep("^(damag|dent|ding)$", dimnames(terms_stem_mtrx)[[2]])], glbObsTrn[, glb_rsp_var], use="pairwise.complete.obs")
    print(subset(glb_post_stem_words_terms_df_lst[[txt_var]], term %in% syn_lst$syns))
    print(subset(get_corpus_terms(tm_map(glb_txt_corpus_lst[[txt_var]], mycombineSynonyms, list(syn_lst), lazy=FALSE)), term == syn_lst$word))
#     cor(terms_stop_mtrx[glbObsAll$.src == "Train", grep("^(damage|dent|ding)$", dimnames(terms_stop_mtrx)[[2]])], glbObsTrn[, glb_rsp_var], use="pairwise.complete.obs")
#     cor(rowSums(terms_stop_mtrx[glbObsAll$.src == "Train", grep("^(damage|dent|ding)$", dimnames(terms_stop_mtrx)[[2]])]), glbObsTrn[, glb_rsp_var], use="pairwise.complete.obs")
}
#chk_comb_cor(syn_lst=list(word="cabl",  syns=c("cabl", "cord")))
#chk_comb_cor(syn_lst=list(word="damag",  syns=c("damag", "dent", "ding")))
#chk_comb_cor(syn_lst=list(word="dent",  syns=c("dent", "ding")))
#chk_comb_cor(syn_lst=list(word="use",  syns=c("use", "usag")))

glb_txt_synonyms <- list()
#glb_txt_synonyms[["<txt_var>"]] <- list(NULL
#     , list(word="<stem1>",  syns=c("<stem1>", "<stem1_2>"))
#                                       )

# options include: "weightTf", "myweightTflog1p", "myweightTfsqrt", "weightTfIdf", "weightBM25"
glb_txt_terms_control <- list(weighting = "weightTfIdf" # : default
                # termFreq selection criteria across obs: tm default: list(global=c(1, Inf))
                    , bounds = list(global = c(1, Inf)) 
                # wordLengths selection criteria: tm default: c(3, Inf)
                    , wordLengths = c(1, Inf) 
                              ) 

glb_txt_cor_var <- glb_rsp_var # : default # or c(<feat>)

# select one from c("union.top.val.cor", "top.cor", "top.val", default: "top.chisq", "sparse")
glbFeatsTextFilter <- "top.chisq" 
glbFeatsTextTermsMax <- rep(10, length(glbFeatsText)) # :default
names(glbFeatsTextTermsMax) <- glbFeatsText

# Text Processing Step: extractAssoc
glbFeatsTextAssocCor <- rep(1, length(glbFeatsText)) # :default 
names(glbFeatsTextAssocCor) <- glbFeatsText

# Remember to use stemmed terms
glb_important_terms <- list()

# Text Processing Step: extractPatterns (ngrams)
glbFeatsTextPatterns <- list()
#glbFeatsTextPatterns[[<txt_var>>]] <- list()
#glbFeatsTextPatterns[[<txt_var>>]] <- c(metropolitan.diary.colon = "Metropolitan Diary:")

# Have to set it even if it is not used
# Properties:
#   numrows(glb_feats_df) << numrows(glbObsFit
#   Select terms that appear in at least 0.2 * O(FP/FN(glbObsOOB)) ???
#       numrows(glbObsOOB) = 1.1 * numrows(glbObsNew) ???
glb_sprs_thresholds <- NULL # or c(<txt_var1> = 0.988, <txt_var2> = 0.970, <txt_var3> = 0.970)

glbFctrMaxUniqVals <- 21 # default: 20
glb_impute_na_data <- TRUE # FALSE # or TRUE
glb_mice_complete.seed <- 144 # or any integer

glb_cluster <- FALSE # : default or TRUE
glb_cluster.seed <- 189 # or any integer
glb_cluster_entropy_var <- glb_rsp_var # c(glb_rsp_var, as.factor(cut(glb_rsp_var, 3)), default: NULL)
glbFeatsTextClusterVarsExclude <- FALSE # default FALSE

glb_interaction_only_feats <- NULL # : default or c(<parent_feat> = "<child_feat>")

glb_nzv_freqCut <- 19 # 19 : caret default
glb_nzv_uniqueCut <- 10 # 10 : caret default

glbRFESizes <- list()
#glbRFESizes[["mdlFamily"]] <- c(4, 8, 16, 32, 64, 67, 68, 69) # Accuracy@69/70 = 0.8258

glbObsFitOutliers <- list()
# If outliers.n >= 10; consider concatenation of interaction vars
# glbObsFitOutliers[["<mdlFamily>"]] <- c(NULL
    # is.na(.rstudent)
    # is.na(.dffits)
    # .hatvalues >= 0.99        
    # -38,167,642 < minmax(.rstudent) < 49,649,823    
#     , <comma-separated-<glb_id_var>>
#                                     )
glbObsTrnOutliers <- list()

# influence.measures: car::outlier; rstudent; dffits; hatvalues; dfbeta; dfbetas
#mdlId <- "RFE.X.glm"; obs_df <- fitobs_df
#mdlId <- "Final.glm"; obs_df <- trnobs_df
#mdlId <- "CSM2.X.glm"; obs_df <- fitobs_df
#print(outliers <- car::outlierTest(glb_models_lst[[mdlId]]$finalModel))
#mdlIdFamily <- paste0(head(unlist(str_split(mdlId, "\\.")), -1), collapse="."); obs_df <- dplyr::filter_(obs_df, interp(~(!(var %in% glbObsFitOutliers[[mdlIdFamily]])), var = as.name(glb_id_var))); model_diags_df <- cbind(obs_df, data.frame(.rstudent=stats::rstudent(glb_models_lst[[mdlId]]$finalModel)), data.frame(.dffits=stats::dffits(glb_models_lst[[mdlId]]$finalModel)), data.frame(.hatvalues=stats::hatvalues(glb_models_lst[[mdlId]]$finalModel)));print(summary(model_diags_df[, c(".rstudent",".dffits",".hatvalues")])); table(cut(model_diags_df$.hatvalues, breaks=c(0.00, 0.98, 0.99, 1.00)))

#print(subset(model_diags_df, is.na(.rstudent))[, glb_id_var])
#print(subset(model_diags_df, is.na(.dffits))[, glb_id_var])
#print(model_diags_df[which.min(model_diags_df$.dffits), ])
#print(subset(model_diags_df, .hatvalues > 0.99)[, glb_id_var])
#dffits_df <- merge(dffits_df, outliers_df, by="row.names", all.x=TRUE); row.names(dffits_df) <- dffits_df$Row.names; dffits_df <- subset(dffits_df, select=-Row.names)
#dffits_df <- merge(dffits_df, glbObsFit, by="row.names", all.x=TRUE); row.names(dffits_df) <- dffits_df$Row.names; dffits_df <- subset(dffits_df, select=-Row.names)
#subset(dffits_df, !is.na(.Bonf.p))

#mdlId <- "CSM.X.glm"; vars <- myextract_actual_feats(row.names(orderBy(reformulate(c("-", paste0(mdlId, ".imp"))), myget_feats_imp(glb_models_lst[[mdlId]])))); 
#model_diags_df <- glb_get_predictions(model_diags_df, mdlId, glb_rsp_var)
#obs_ix <- row.names(model_diags_df) %in% names(outliers$rstudent)[1]
#obs_ix <- which(is.na(model_diags_df$.rstudent))
#obs_ix <- which(is.na(model_diags_df$.dffits))
#myplot_parcoord(obs_df=model_diags_df[, c(glb_id_var, glbFeatsCategory, ".rstudent", ".dffits", ".hatvalues", glb_rsp_var, paste0(glb_rsp_var, mdlId), vars[1:min(20, length(vars))])], obs_ix=obs_ix, id_var=glb_id_var, category_var=glbFeatsCategory)

#model_diags_df[row.names(model_diags_df) %in% names(outliers$rstudent)[c(1:2)], ]
#ctgry_diags_df <- model_diags_df[model_diags_df[, glbFeatsCategory] %in% c("Unknown#0"), ]
#myplot_parcoord(obs_df=ctgry_diags_df[, c(glb_id_var, glbFeatsCategory, ".rstudent", ".dffits", ".hatvalues", glb_rsp_var, "startprice.log10.predict.RFE.X.glmnet", indep_vars[1:20])], obs_ix=row.names(ctgry_diags_df) %in% names(outliers$rstudent)[1], id_var=glb_id_var, category_var=glbFeatsCategory)
#table(glbObsFit[model_diags_df[, glbFeatsCategory] %in% c("iPad1#1"), "startprice.log10.cut.fctr"])
#glbObsFit[model_diags_df[, glbFeatsCategory] %in% c("iPad1#1"), c(glb_id_var, "startprice")]

# No outliers & .dffits == NaN
#myplot_parcoord(obs_df=model_diags_df[, c(glb_id_var, glbFeatsCategory, glb_rsp_var, "startprice.log10.predict.RFE.X.glmnet", indep_vars[1:10])], obs_ix=seq(1:nrow(model_diags_df))[is.na(model_diags_df$.dffits)], id_var=glb_id_var, category_var=glbFeatsCategory)

# Modify mdlId to (build & extract) "<FamilyId>#<Fit|Trn>#<caretMethod>#<preProc1.preProc2>#<samplingMethod>"
glb_models_lst <- list(); glb_models_df <- data.frame()
# Regression
if (glb_is_regression) {
    glbMdlMethods <- c(NULL
        # deterministic
            #, "lm", # same as glm
            , "glm", "bayesglm", "glmnet"
            , "rpart"
        # non-deterministic
            , "gbm", "rf" 
        # Unknown
            , "nnet" , "avNNet" # runs 25 models per cv sample for tunelength=5
            , "svmLinear", "svmLinear2"
            , "svmPoly" # runs 75 models per cv sample for tunelength=5
            , "svmRadial" 
            , "earth"
            , "bagEarth" # Takes a long time
        )
} else
# Classification - Add ada (auto feature selection)
    if (glb_is_binomial)
        glbMdlMethods <- c(NULL
        # deterministic                     
            , "bagEarth" # Takes a long time        
            , "glm", "bayesglm", "glmnet"
            , "nnet"
            , "rpart"
        # non-deterministic        
            , "gbm"
            , "avNNet" # runs 25 models per cv sample for tunelength=5      
            , "rf"
        # Unknown
            , "lda", "lda2"
                # svm models crash when predict is called -> internal to kernlab it should call predict without .outcome
            , "svmLinear", "svmLinear2"
            , "svmPoly" # runs 75 models per cv sample for tunelength=5
            , "svmRadial" 
            , "earth"
        ) else
        glbMdlMethods <- c(NULL
        # non-deterministic 
            , "rf"       
        # Unknown
            , "gbm", "rpart"
        )

glb_mdl_family_lst <- list(); glb_mdl_feats_lst <- list()
# family: Choose from c("RFE.X", "CSM.X", "All.X", "Best.Interact")
#   methods: Choose from c(NULL, <method>, glbMdlMethods) 
#glb_mdl_family_lst[["RFE.X"]] <- c("glmnet", "glm") # non-NULL list is mandatory
glb_mdl_family_lst[["All.X"]] <- "glmnet" # non-NULL list is mandatory
#glb_mdl_family_lst[["Best.Interact"]] <- "glmnet" # non-NULL list is mandatory

# Check if interaction features make RFE better
# glb_mdl_family_lst[["CSM.X"]] <- setdiff(glbMdlMethods, c("lda", "lda2")) # crashing due to category:.clusterid ??? #c("glmnet", "glm") # non-NULL list is mandatory
# glb_mdl_feats_lst[["CSM.X"]] <- c(NULL
#     , <comma-separated-features-vector>
#                                   )
# dAFeats.CSM.X %<d-% c(NULL
#     # Interaction feats up to varImp(RFE.X.glmnet) >= 50
#     , <comma-separated-features-vector>
#     , setdiff(myextract_actual_feats(predictors(rfe_fit_results)), c(NULL
#                , <comma-separated-features-vector>
#                                                                       ))    
#                                   )
# glb_mdl_feats_lst[["CSM.X"]] <- "%<d-% dAFeats.CSM.X"

# Check if tuning parameters make fit better; make it mdlFamily customizable ?
glbMdlTuneParams <- data.frame()
# When glmnet crashes at model$grid with error: ???
glmnetTuneParams <- rbind(data.frame()
                        ,data.frame(parameter = "alpha",  vals = "0.100 0.325 0.550 0.775 1.000")
                        ,data.frame(parameter = "lambda", vals = "9.342e-02")    
                        )
glbMdlTuneParams <- myrbind_df(glbMdlTuneParams,
                               cbind(data.frame(mdlId = "Max.cor.Y.Time.Lag##rcv#glmnet"),
                                     glmnetTuneParams))

    #avNNet    
    #   size=[1] 3 5 7 9; decay=[0] 1e-04 0.001  0.01   0.1; bag=[FALSE]; RMSE=1.3300906 

    #bagEarth
    #   degree=1 [2] 3; nprune=64 128 256 512 [1024]; RMSE=0.6486663 (up)
# glbMdlTuneParams <- myrbind_df(glbMdlTuneParams, rbind(data.frame()
#     ,data.frame(method = "bagEarth", parameter = "nprune", vals = "256")
#     ,data.frame(method = "bagEarth", parameter = "degree", vals = "2")    
# ))

    #earth 
    #   degree=[1]; nprune=2  [9] 17 25 33; RMSE=0.1334478
    
    #gbm 
    #   shrinkage=0.05 [0.10] 0.15 0.20 0.25; n.trees=100 150 200 [250] 300; interaction.depth=[1] 2 3 4 5; n.minobsinnode=[10]; RMSE=0.2008313     
# glbMdlTuneParams <- myrbind_df(glbMdlTuneParams, rbind(data.frame()
#     ,data.frame(method = "gbm", parameter = "shrinkage", min = 0.05, max = 0.25, by = 0.05)
#     ,data.frame(method = "gbm", parameter = "n.trees", min = 100, max = 300, by = 50)
#     ,data.frame(method = "gbm", parameter = "interaction.depth", min = 1, max = 5, by = 1)
#     ,data.frame(method = "gbm", parameter = "n.minobsinnode", min = 10, max = 10, by = 10)
#     #seq(from=0.05,  to=0.25, by=0.05)
# ))

    #glmnet
    #   alpha=0.100 [0.325] 0.550 0.775 1.000; lambda=0.0005232693 0.0024288010 0.0112734954 [0.0523269304] 0.2428800957; RMSE=0.6164891
# glbMdlTuneParams <- myrbind_df(glbMdlTuneParams, rbind(data.frame()
#     ,data.frame(method = "glmnet", parameter = "alpha", vals = "0.550 0.775 0.8875 0.94375 1.000")
#     ,data.frame(method = "glmnet", parameter = "lambda", vals = "9.858855e-05 0.0001971771 0.0009152152 0.0042480525 0.0197177130")    
# ))

    #nnet    
    #   size=3 5 [7] 9 11; decay=0.0001 0.001 0.01 [0.1] 0.2; RMSE=0.9287422
# glbMdlTuneParams <- myrbind_df(glbMdlTuneParams, rbind(data.frame()
#     ,data.frame(method = "nnet", parameter = "size", vals = "3 5 7 9 11")
#     ,data.frame(method = "nnet", parameter = "decay", vals = "0.0001 0.0010 0.0100 0.1000 0.2000")    
# ))

    #rf # Don't bother; results are not deterministic
    #       mtry=2  35  68 [101] 134; RMSE=0.1339974
# glbMdlTuneParams <- myrbind_df(glbMdlTuneParams, rbind(data.frame()
#     ,data.frame(method = "rf", parameter = "mtry", vals = "2 5 9 13 17")
# ))

    #rpart 
    #   cp=0.020 [0.025] 0.030 0.035 0.040; RMSE=0.1770237
# glbMdlTuneParams <- myrbind_df(glbMdlTuneParams, rbind(data.frame()    
#     ,data.frame(method = "rpart", parameter = "cp", vals = "0.004347826 0.008695652 0.017391304 0.021739130 0.034782609")
# ))
    
    #svmLinear
    #   C=0.01 0.05 [0.10] 0.50 1.00 2.00 3.00 4.00; RMSE=0.1271318; 0.1296718
# glbMdlTuneParams <- myrbind_df(glbMdlTuneParams, rbind(data.frame()
#     ,data.frame(method = "svmLinear", parameter = "C", vals = "0.01 0.05 0.1 0.5 1")
# ))

    #svmLinear2    
    #   cost=0.0625 0.1250 [0.25] 0.50 1.00; RMSE=0.1276354 
# glbMdlTuneParams <- myrbind_df(glbMdlTuneParams, rbind(data.frame()
#     ,data.frame(method = "svmLinear2", parameter = "cost", vals = "0.0625 0.125 0.25 0.5 1")
# ))

    #svmPoly    
    #   degree=[1] 2 3 4 5; scale=0.01 0.05 [0.1] 0.5 1; C=0.50 1.00 [2.00] 3.00 4.00; RMSE=0.1276130
# glbMdlTuneParams <- myrbind_df(glbMdlTuneParams, rbind(data.frame()
#     ,data.frame(method="svmPoly", parameter="degree", min=1, max=5, by=1) #seq(1, 5, 1)
#     ,data.frame(method="svmPoly", parameter="scale", vals="0.01, 0.05, 0.1, 0.5, 1")
#     ,data.frame(method="svmPoly", parameter="C", vals="0.50, 1.00, 2.00, 3.00, 4.00")    
# ))

    #svmRadial
    #   sigma=[0.08674323]; C=0.25 0.50 1.00 [2.00] 4.00; RMSE=0.1614957
    
#glb2Sav(); all.equal(sav_models_df, glb_models_df)
    
glb_preproc_methods <- NULL
#     c("YeoJohnson", "center.scale", "range", "pca", "ica", "spatialSign")

# Baseline prediction model feature(s)
glb_Baseline_mdl_var <- NULL # or c("<feat>")

glbMdlMetric_terms <- NULL # or matrix(c(
#                               0,1,2,3,4,
#                               2,0,1,2,3,
#                               4,2,0,1,2,
#                               6,4,2,0,1,
#                               8,6,4,2,0
#                           ), byrow=TRUE, nrow=5)
glbMdlMetricSummary <- NULL # or "<metric_name>"
glbMdlMetricMaximize <- NULL # or FALSE (TRUE is not the default for both classification & regression) 
glbMdlMetricSummaryFn <- NULL # or function(data, lev=NULL, model=NULL) {
#     confusion_mtrx <- t(as.matrix(confusionMatrix(data$pred, data$obs)))
#     #print(confusion_mtrx)
#     #print(confusion_mtrx * glbMdlMetric_terms)
#     metric <- sum(confusion_mtrx * glbMdlMetric_terms) / nrow(data)
#     names(metric) <- glbMdlMetricSummary
#     return(metric)
# }

glbMdlCheckRcv <- FALSE # Turn it on when needed; otherwise takes long time
glb_rcv_n_folds <- 3 # or NULL
glb_rcv_n_repeats <- 3 # or NULL

glb_clf_proba_threshold <- NULL # 0.5

# Model selection criteria
if (glb_is_regression)
    glbMdlMetricsEval <- c("min.RMSE.OOB", "max.R.sq.OOB", "max.Adj.R.sq.fit", "min.RMSE.fit")
    #glbMdlMetricsEval <- c("min.RMSE.fit", "max.R.sq.fit", "max.Adj.R.sq.fit")    
if (glb_is_classification) {
    if (glb_is_binomial)
        glbMdlMetricsEval <- 
            c("max.Accuracy.OOB", "max.AUCROCR.OOB", "max.AUCpROC.OOB", "min.aic.fit", "max.Accuracy.fit") else        
        glbMdlMetricsEval <- c("max.Accuracy.OOB", "max.Kappa.OOB")
}

# select from NULL [no ensemble models], "auto" [all models better than MFO or Baseline], c(mdl_ids in glb_models_lst) [Typically top-rated models in auto]
glb_mdl_ensemble <- NULL
#     "%<d-% setdiff(mygetEnsembleAutoMdlIds(), 'CSM.X.rf')" 
#     c(<comma-separated-mdlIds>
#      )

# Only for classifications; for regressions remove "(.*)\\.prob" form the regex
# tmp_fitobs_df <- glbObsFit[, grep(paste0("^", gsub(".", "\\.", mygetPredictIds$value, fixed = TRUE), "CSM\\.X\\.(.*)\\.prob"), names(glbObsFit), value = TRUE)]; cor_mtrx <- cor(tmp_fitobs_df); cor_vctr <- sort(cor_mtrx[row.names(orderBy(~-Overall, varImp(glb_models_lst[["Ensemble.repeatedcv.glmnet"]])$imp))[1], ]); summary(cor_vctr); cor_vctr
#ntv.glm <- glm(reformulate(indep_vars, glb_rsp_var), family = "binomial", data = glbObsFit)
#step.glm <- step(ntv.glm)

glb_sel_mdl_id <- "All.X##rcv#glmnet" #select from c(NULL, "All.X##rcv#glmnet", "RFE.X##rcv#glmnet", <mdlId>)
glb_fin_mdl_id <- NULL #select from c(NULL, glb_sel_mdl_id)

glb_dsp_cols <- c(glb_id_var, glbFeatsCategory, glb_rsp_var
#               List critical cols excl. glb_id_var, glbFeatsCategory & glb_rsp_var
                  )

# Output specs
glbOutDataVizFname <- "NYTBlogs3_obsall.csv" # choose from c(NULL, "NYTBlogs3_obsall.csv")
glb_out_obs <- NULL # select from c(NULL : default to "new", "all", "new", "trn")
glb_out_vars_lst <- list()
# glb_id_var will be the first output column, by default
glb_out_vars_lst[["Probability1"]] <- 
    "%<d-% mygetPredictIds(glb_rsp_var, glb_fin_mdl_id)$prob"
# glb_out_vars_lst[[glb_rsp_var_raw]] <- glb_rsp_var_raw
# glb_out_vars_lst[[paste0(head(unlist(strsplit(mygetPredictIds$value, "")), -1), collapse = "")]] <-

glbOutStackFnames <- NULL #: default
    # c("ebayipads_txt_assoc1_out_bid1_stack.csv") # manual stack
    # c("ebayipads_finmdl_bid1_out_nnet_1.csv") # universal stack
glb_out_pfx <- "NYTBlogs3_feat_PubDate_"
glb_save_envir <- FALSE # or TRUE

# Depict process
glb_analytics_pn <- petrinet(name = "glb_analytics_pn",
                        trans_df = data.frame(id = 1:6,
    name = c("data.training.all","data.new",
           "model.selected","model.final",
           "data.training.all.prediction","data.new.prediction"),
    x=c(   -5,-5,-15,-25,-25,-35),
    y=c(   -5, 5,  0,  0, -5,  5)
                        ),
                        places_df=data.frame(id=1:4,
    name=c("bgn","fit.data.training.all","predict.data.new","end"),
    x=c(   -0,   -20,                    -30,               -40),
    y=c(    0,     0,                      0,                 0),
    M0=c(   3,     0,                      0,                 0)
                        ),
                        arcs_df=data.frame(
    begin=c("bgn","bgn","bgn",        
            "data.training.all","model.selected","fit.data.training.all",
            "fit.data.training.all","model.final",    
            "data.new","predict.data.new",
            "data.training.all.prediction","data.new.prediction"),
    end  =c("data.training.all","data.new","model.selected",
            "fit.data.training.all","fit.data.training.all","model.final",
            "data.training.all.prediction","predict.data.new",
            "predict.data.new","data.new.prediction",
            "end","end")
                        ))
#print(ggplot.petrinet(glb_analytics_pn))
print(ggplot.petrinet(glb_analytics_pn) + coord_flip())
```

```
## Loading required package: grid
```

![](NYTBlogs3_feat_PubDate_files/figure-html/set_global_options-1.png) 

```r
glb_analytics_avl_objs <- NULL

glb_chunks_df <- myadd_chunk(NULL, "import.data")
```

```
##         label step_major step_minor label_minor    bgn end elapsed
## 1 import.data          1          0           0 11.375  NA      NA
```

## Step `1.0: import data`
#### chunk option: eval=<r condition>

```
## [1] "Reading file ./data/NYTimesBlogTrain.csv..."
## [1] "dimensions of data in ./data/NYTimesBlogTrain.csv: 6,532 rows x 10 cols"
##   NewsDesk      SectionName SubsectionName
## 1 Business Crosswords/Games               
## 2  Culture             Arts               
## 3 Business     Business Day       Dealbook
## 4 Business     Business Day       Dealbook
## 5  Science           Health               
## 6  Science           Health               
##                                            Headline
## 1                                  More School Daze
## 2      New 96-Page Murakami Work Coming in December
## 3 Public Pension Funds Stay Mum on Corporate Expats
## 4                             Boot Camp for Bankers
## 5                     Of Little Help to Older Knees
## 6                     A Benefit of Legal Marijuana 
##                                                                                                                                                                                                                           Snippet
## 1                                                                                                                                                                  A puzzle from Ethan Cooper that reminds me that a bill is due.
## 2                                                                            The Strange Library will arrive just three and a half months after Mr. Murakamis latest novel, Colorless Tsukuru Tazaki and His Years of Pilgrimage.
## 3                                      Public pension funds have major stakes in American companies moving overseas to cut their tax bills. But they are saying little about the strategy, which could hurt the nations tax base.
## 4                                                                         As they struggle to find new business to bolster sluggish earnings, banks consider the nations 25 million veterans and service members ideal customers.
## 5                                         Middle-aged and older patients are unlikely to benefit in the long term from surgery to repair tears in the meniscus, pads of cartilage in the knee, a new review of studies has found.
## 6 A new study has found evidence that legal access to marijuana is associated with fewer opioid overdose deaths, but researchers said their findings should not be used as the basis for the wide adoption of legalized cannabis.
##                                                                                                                                                                                                                          Abstract
## 1                                                                                                                                                                  A puzzle from Ethan Cooper that reminds me that a bill is due.
## 2                                                                            The Strange Library will arrive just three and a half months after Mr. Murakamis latest novel, Colorless Tsukuru Tazaki and His Years of Pilgrimage.
## 3                                      Public pension funds have major stakes in American companies moving overseas to cut their tax bills. But they are saying little about the strategy, which could hurt the nations tax base.
## 4                                                                         As they struggle to find new business to bolster sluggish earnings, banks consider the nations 25 million veterans and service members ideal customers.
## 5                                         Middle-aged and older patients are unlikely to benefit in the long term from surgery to repair tears in the meniscus, pads of cartilage in the knee, a new review of studies has found.
## 6 A new study has found evidence that legal access to marijuana is associated with fewer opioid overdose deaths, but researchers said their findings should not be used as the basis for the wide adoption of legalized cannabis.
##   WordCount             PubDate Popular UniqueID
## 1       508 2014-09-01 22:00:09       1        1
## 2       285 2014-09-01 21:14:07       0        2
## 3      1211 2014-09-01 21:05:36       0        3
## 4      1405 2014-09-01 20:43:34       1        4
## 5       181 2014-09-01 18:58:51       1        5
## 6       245 2014-09-01 18:52:22       1        6
##      NewsDesk      SectionName SubsectionName
## 226    Styles                                
## 995                                          
## 3327                                         
## 4753                Multimedia               
## 4802 Business Crosswords/Games               
## 6463   TStyle                                
##                                                   Headline
## 226  For Tavi Gevinson, Fashion Takes a Back Seat, for Now
## 995          Reconsidering What to Call an Extremist Group
## 3327     Clinton's Diagnosis of What's Wrong With Politics
## 4753       'Off Color' and on Target About Race in America
## 4802                      Daniel Finkel's Circle-Toss Game
## 6463                                     Entering the Void
##                                                                                                                                                                            Snippet
## 226                                                Tavi Gevinson, the teenage fashion star turned Broadway actress, wont be much of a player at New York Fashion Week this season.
## 995                                                    Editors have decided to adjust how The Times refer to an Islamic extremist group that controls territory in Syria and Iraq.
## 3327 Hillary Rodham Clinton continued to laugh off questions about her presidential aspirations on Tuesday, but she did shed some light on what she thinks is wrong in Washington.
## 4753              Off Color, a New York Times video series, looks at how artists of color are making sharp social commentary about race in America through comedy and performance.
## 4802                                                                                            By math educator Daniel Finkel, a puzzle thats childs play. Can you figure it out?
## 6463                      The Spanish artist Miquel Barcel closely examines the basic materials of life in response to Edward Hirsch questioning his own belief in a higher power.
##                                                                                                                                                                           Abstract
## 226                                                Tavi Gevinson, the teenage fashion star turned Broadway actress, wont be much of a player at New York Fashion Week this season.
## 995                                                    Editors have decided to adjust how The Times refer to an Islamic extremist group that controls territory in Syria and Iraq.
## 3327 Hillary Rodham Clinton continued to laugh off questions about her presidential aspirations on Tuesday, but she did shed some light on what she thinks is wrong in Washington.
## 4753              Off Color, a New York Times video series, looks at how artists of color are making sharp social commentary about race in America through comedy and performance.
## 4802                                                                                            By math educator Daniel Finkel, a puzzle thats childs play. Can you figure it out?
## 6463                      The Spanish artist Miquel Barcel closely examines the basic materials of life in response to Edward Hirsch questioning his own belief in a higher power.
##      WordCount             PubDate Popular UniqueID
## 226        459 2014-09-04 16:55:57       0      226
## 995        301 2014-09-15 16:05:13       0      995
## 3327       236 2014-10-14 14:45:51       0     3327
## 4753       393 2014-11-02 05:00:13       0     4753
## 4802      1628 2014-11-03 12:00:04       1     4802
## 6463       264 2014-11-27 12:00:09       0     6463
##      NewsDesk SectionName  SubsectionName
## 6527  Foreign                            
## 6528              Opinion Room For Debate
## 6529  Foreign                            
## 6530   TStyle                            
## 6531           Multimedia                
## 6532 Business                            
##                                                                        Headline
## 6527                                     1914: Russians Dominate in East Poland
## 6528                                             Finding a Secretary of Defense
## 6529                         1889: Metropolitan Opera House Reopens in New York
## 6530                         The Daily Gift: Picasso Plates for Creative Dining
## 6531                                          Racing From New York to Barcelona
## 6532 Math Anxiety: Why Hollywood Makes Robots of Alan Turing and Other Geniuses
##                                                                                                                                                                                             Snippet
## 6527                                                                                                      From the International Herald Tribune archives: Russians dominate in East Poland in 1914.
## 6528                                                                                             If Chuck Hagel isn't the right Pentagon chief to respond to an onslaught of global crises, who is?
## 6529                                                                                      From the International Herald Tribune archives: The Metropolitan Opera House reopens in New York in 1889.
## 6530                                                                                                                      Each day until Christmas, the editors of T share a new holiday gift idea.
## 6531                                                      A sailboat race from New York to Barcelona was the setting for a thrilling  and sometimes terrifying  video about this challenging sport.
## 6532 The visionary who stares at formulas written on walls or mirrors  or better yet, thin air  has become a Hollywood trope. So has the depiction of the genius who cant connect with real people.
##                                                                                                                                                                                            Abstract
## 6527                                                                                                      From the International Herald Tribune archives: Russians dominate in East Poland in 1914.
## 6528                                                                                             If Chuck Hagel isn't the right Pentagon chief to respond to an onslaught of global crises, who is?
## 6529                                                                                      From the International Herald Tribune archives: The Metropolitan Opera House reopens in New York in 1889.
## 6530                                                                                                                      Each day until Christmas, the editors of T share a new holiday gift idea.
## 6531                                                      A sailboat race from New York to Barcelona was the setting for a thrilling  and sometimes terrifying  video about this challenging sport.
## 6532 The visionary who stares at formulas written on walls or mirrors  or better yet, thin air  has become a Hollywood trope. So has the depiction of the genius who cant connect with real people.
##      WordCount             PubDate Popular UniqueID
## 6527       176 2014-11-30 13:48:40       0     6527
## 6528      1597 2014-11-30 13:27:23       0     6528
## 6529       214 2014-11-30 09:44:57       0     6529
## 6530        61 2014-11-30 09:00:43       0     6530
## 6531       441 2014-11-30 09:00:22       0     6531
## 6532       921 2014-11-30 07:00:40       0     6532
## 'data.frame':	6532 obs. of  10 variables:
##  $ NewsDesk      : chr  "Business" "Culture" "Business" "Business" ...
##  $ SectionName   : chr  "Crosswords/Games" "Arts" "Business Day" "Business Day" ...
##  $ SubsectionName: chr  "" "" "Dealbook" "Dealbook" ...
##  $ Headline      : chr  "More School Daze" "New 96-Page Murakami Work Coming in December" "Public Pension Funds Stay Mum on Corporate Expats" "Boot Camp for Bankers" ...
##  $ Snippet       : chr  "A puzzle from Ethan Cooper that reminds me that a bill is due." "The Strange Library will arrive just three and a half months after Mr. Murakamis latest novel, Colorless Tsukuru Tazaki and His"| __truncated__ "Public pension funds have major stakes in American companies moving overseas to cut their tax bills. But they are saying little"| __truncated__ "As they struggle to find new business to bolster sluggish earnings, banks consider the nations 25 million veterans and service "| __truncated__ ...
##  $ Abstract      : chr  "A puzzle from Ethan Cooper that reminds me that a bill is due." "The Strange Library will arrive just three and a half months after Mr. Murakamis latest novel, Colorless Tsukuru Tazaki and His"| __truncated__ "Public pension funds have major stakes in American companies moving overseas to cut their tax bills. But they are saying little"| __truncated__ "As they struggle to find new business to bolster sluggish earnings, banks consider the nations 25 million veterans and service "| __truncated__ ...
##  $ WordCount     : int  508 285 1211 1405 181 245 258 893 1077 188 ...
##  $ PubDate       : chr  "2014-09-01 22:00:09" "2014-09-01 21:14:07" "2014-09-01 21:05:36" "2014-09-01 20:43:34" ...
##  $ Popular       : int  1 0 0 1 1 1 0 1 1 0 ...
##  $ UniqueID      : int  1 2 3 4 5 6 7 8 9 10 ...
##  - attr(*, "comment")= chr "glbObsTrn"
## NULL
```

```
## [1] "Reading file ./data/NYTimesBlogTest.csv..."
## [1] "dimensions of data in ./data/NYTimesBlogTest.csv: 1,870 rows x 9 cols"
##   NewsDesk      SectionName SubsectionName
## 1  Culture                                
## 2  Culture             Arts               
## 3 Business Crosswords/Games               
## 4 Business     Business Day       Dealbook
## 5  Science           Health               
## 6  Science           Health               
##                                                             Headline
## 1                                         'Birdman' Tops the Gothams
## 2                     'Sleepy Hollow' Recap: A Not-So-Shocking Death
## 3                                        Drinking Buddy For Falstaff
## 4 Encouraging Public Service, Through Wall Street's 'Revolving Door'
## 5                           Therapy Prevents Repeat Suicide Attempts
## 6                                            Hoping for a Good Death
##                                                                                                                                                 Snippet
## 1                                                    The backstage tale won two awards; Citizenfour, the Edward Snowden documentary, was also a winner.
## 2                                                                      In the fall season finale, a question of where the series has many places to go.
## 3                                                                                                       In which Timothy Polin reveals his potty mouth.
## 4 The debate about pay for Wall Street executives who take government jobs appears to be based more on a populist shakedown than on good public policy.
## 5                                                                Short-term psychotherapy may be an effective way to prevent repeated suicide attempts.
## 6                          What I hadnt considered before my fathers heart attack was the precise meaning of not wanting to live hooked up to machines.
##                                                                                                                                                Abstract
## 1                                                    The backstage tale won two awards; Citizenfour, the Edward Snowden documentary, was also a winner.
## 2                                                                      In the fall season finale, a question of where the series has many places to go.
## 3                                                                                                       In which Timothy Polin reveals his potty mouth.
## 4 The debate about pay for Wall Street executives who take government jobs appears to be based more on a populist shakedown than on good public policy.
## 5                                                                Short-term psychotherapy may be an effective way to prevent repeated suicide attempts.
## 6                          What I hadnt considered before my fathers heart attack was the precise meaning of not wanting to live hooked up to machines.
##   WordCount             PubDate UniqueID
## 1       111 2014-12-01 22:45:24     6533
## 2       558 2014-12-01 22:01:34     6534
## 3       788 2014-12-01 22:00:26     6535
## 4       915 2014-12-01 21:04:13     6536
## 5       213 2014-12-01 19:13:20     6537
## 6       938 2014-12-01 19:05:12     6538
##     NewsDesk      SectionName SubsectionName
## 3   Business Crosswords/Games               
## 334     OpEd          Opinion               
## 725   TStyle                                
## 732 Business     Business Day       Dealbook
## 752 Business     Business Day       Dealbook
## 864                                         
##                                                            Headline
## 3                                       Drinking Buddy For Falstaff
## 334 Facts & Figures: America&rsquo;s Unique Take on Maternity Leave
## 725                               Ansel Elgort Buttons Up in Brioni
## 732      A Shake-Up as the Financial World Infiltrates Philanthropy
## 752    Coupang, a South Korean E-Commerce Site, Raises $300 Million
## 864                                               Today in Politics
##                                                                                                                                                 Snippet
## 3                                                                                                       In which Timothy Polin reveals his potty mouth.
## 334                                                                                In the U.S., paid parental leave is more of a perk than a guarantee.
## 725                                                        The actor brought a tinge of youthfulness to the classic Italian houses retro-tailored look.
## 732 Donor-advised funds help investors get deductions for charitable donations in one year, but society doesnt get the benefit of the money right away.
## 752                                 The latest financing round underscores Coupangs maturity and its ambitions to one day be a publicly traded company.
## 864          The 113th Congress is concluding with partisan brinksmanship and one last mad scramble for votes to pass a $1.1 trillion spending package.
##                                                                                                                                                Abstract
## 3                                                                                                       In which Timothy Polin reveals his potty mouth.
## 334                                                                                In the U.S., paid parental leave is more of a perk than a guarantee.
## 725                                                        The actor brought a tinge of youthfulness to the classic Italian houses retro-tailored look.
## 732 Donor-advised funds help investors get deductions for charitable donations in one year, but society doesnt get the benefit of the money right away.
## 752                                 The latest financing round underscores Coupangs maturity and its ambitions to one day be a publicly traded company.
## 864          The 113th Congress is concluding with partisan brinksmanship and one last mad scramble for votes to pass a $1.1 trillion spending package.
##     WordCount             PubDate UniqueID
## 3         788 2014-12-01 22:00:26     6535
## 334       160 2014-12-04 11:45:20     6866
## 725        89 2014-12-10 12:30:47     7257
## 732      1172 2014-12-10 12:00:38     7264
## 752       353 2014-12-10 08:30:41     7284
## 864      1544 2014-12-11 07:09:25     7396
##      NewsDesk   SectionName SubsectionName
## 1865                                      
## 1866 Business    Technology               
## 1867    Metro N.Y. / Region               
## 1868             Multimedia               
## 1869  Foreign         World   Asia Pacific
## 1870  Science        Health               
##                                                       Headline
## 1865                                         Today in Politics
## 1866                         Uber Suspends Operations in Spain
## 1867                         New York Today: The Year in News 
## 1868                   New Year, Old Memories, in Times Square
## 1869 Hong Kong Police Criticized After 14-Year-Old's Detention
## 1870          The Super-Short Workout and Other Fitness Trends
##                                                                                                                                                                                                                                                   Snippet
## 1865                                                                                                               House Republicans are ending the year on a defensive note over Representative Steve Scalises 2002 speech to a white supremacist group.
## 1866                                                                              In a first in the growing pushback against Ubers global expansion, a judges ruling barred telecommunications operators and banks from supporting the companys services.
## 1867                                                                                                                                                              Wednesday: The most read stories of 2014, teeth-chattering cold, and its New Years Eve.
## 1868                                                                         What happens when you combine Burning Man, Independence Day fireworks, the last day of school and a full-contact Black Friday sale-a-bration? New Years Eve in Times Square.
## 1869 The authorities have been accused of trying to intimidate young pro-democracy protesters and their families after a 14-year-old girl was detained on suspicion of drawing flowers in chalk near government headquarters and sent to a juvenile home.
## 1870                                                                                                                 The big story in exercise science this year was the super-short workout, although many other fitness-related themes emerged in 2014.
##                                                                                                                                                                                                                                                  Abstract
## 1865                                                                                                               House Republicans are ending the year on a defensive note over Representative Steve Scalises 2002 speech to a white supremacist group.
## 1866                                                                              In a first in the growing pushback against Ubers global expansion, a judges ruling barred telecommunications operators and banks from supporting the companys services.
## 1867                                                                                                                                                              Wednesday: The most read stories of 2014, teeth-chattering cold, and its New Years Eve.
## 1868                                                                         What happens when you combine Burning Man, Independence Day fireworks, the last day of school and a full-contact Black Friday sale-a-bration? New Years Eve in Times Square.
## 1869 The authorities have been accused of trying to intimidate young pro-democracy protesters and their families after a 14-year-old girl was detained on suspicion of drawing flowers in chalk near government headquarters and sent to a juvenile home.
## 1870                                                                                                                 The big story in exercise science this year was the super-short workout, although many other fitness-related themes emerged in 2014.
##      WordCount             PubDate UniqueID
## 1865      1616 2014-12-31 07:03:46     8397
## 1866       292 2014-12-31 06:09:32     8398
## 1867      1010 2014-12-31 06:06:58     8399
## 1868       387 2014-12-31 05:00:19     8400
## 1869       717 2014-12-31 04:16:29     8401
## 1870       818 2014-12-31 00:01:10     8402
## 'data.frame':	1870 obs. of  9 variables:
##  $ NewsDesk      : chr  "Culture" "Culture" "Business" "Business" ...
##  $ SectionName   : chr  "" "Arts" "Crosswords/Games" "Business Day" ...
##  $ SubsectionName: chr  "" "" "" "Dealbook" ...
##  $ Headline      : chr  "'Birdman' Tops the Gothams" "'Sleepy Hollow' Recap: A Not-So-Shocking Death" "Drinking Buddy For Falstaff" "Encouraging Public Service, Through Wall Street's 'Revolving Door'" ...
##  $ Snippet       : chr  "The backstage tale won two awards; Citizenfour, the Edward Snowden documentary, was also a winner." "In the fall season finale, a question of where the series has many places to go." "In which Timothy Polin reveals his potty mouth." "The debate about pay for Wall Street executives who take government jobs appears to be based more on a populist shakedown than "| __truncated__ ...
##  $ Abstract      : chr  "The backstage tale won two awards; Citizenfour, the Edward Snowden documentary, was also a winner." "In the fall season finale, a question of where the series has many places to go." "In which Timothy Polin reveals his potty mouth." "The debate about pay for Wall Street executives who take government jobs appears to be based more on a populist shakedown than "| __truncated__ ...
##  $ WordCount     : int  111 558 788 915 213 938 1336 2644 752 99 ...
##  $ PubDate       : chr  "2014-12-01 22:45:24" "2014-12-01 22:01:34" "2014-12-01 22:00:26" "2014-12-01 21:04:13" ...
##  $ UniqueID      : int  6533 6534 6535 6536 6537 6538 6539 6540 6541 6542 ...
##  - attr(*, "comment")= chr "glbObsNew"
## NULL
```

```
## [1] "Partition stats:"
```

```
## Loading required package: sqldf
## Loading required package: gsubfn
## Loading required package: proto
## Loading required package: RSQLite
## Loading required package: DBI
## Loading required package: tcltk
```

```
##   Popular  .src   .n
## 1       0 Train 5439
## 2      NA  Test 1870
## 3       1 Train 1093
##   Popular  .src   .n
## 1       0 Train 5439
## 2      NA  Test 1870
## 3       1 Train 1093
```

![](NYTBlogs3_feat_PubDate_files/figure-html/import.data-1.png) 

```
##    .src   .n
## 1 Train 6532
## 2  Test 1870
```

```
## Loading required package: dplyr
## 
## Attaching package: 'dplyr'
## 
## The following objects are masked from 'package:stats':
## 
##     filter, lag
## 
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
## 
## Loading required package: lazyeval
## Loading required package: gdata
## gdata: read.xls support for 'XLS' (Excel 97-2004) files ENABLED.
## 
## gdata: read.xls support for 'XLSX' (Excel 2007+) files ENABLED.
## 
## Attaching package: 'gdata'
## 
## The following objects are masked from 'package:dplyr':
## 
##     combine, first, last
## 
## The following object is masked from 'package:stats':
## 
##     nobs
## 
## The following object is masked from 'package:utils':
## 
##     object.size
```

```
## [1] "Found 0 duplicates by all features:"
```

```
## NULL
```

```
##          label step_major step_minor label_minor    bgn    end elapsed
## 1  import.data          1          0           0 11.375 21.011   9.636
## 2 inspect.data          2          0           0 21.011     NA      NA
```

## Step `2.0: inspect data`

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

```
## Warning: Removed 1870 rows containing non-finite values (stat_bin).
```

```
## Loading required package: reshape2
```

![](NYTBlogs3_feat_PubDate_files/figure-html/inspect.data-1.png) 

```
##       Popular.0 Popular.1 Popular.NA
## Test         NA        NA       1870
## Train      5439      1093         NA
##       Popular.0 Popular.1 Popular.NA
## Test         NA        NA          1
## Train 0.8326699 0.1673301         NA
## [1] "numeric data missing in glbObsAll: "
## Popular 
##    1870 
## [1] "numeric data w/ 0s in glbObsAll: "
## WordCount   Popular 
##       109      5439 
## [1] "numeric data w/ Infs in glbObsAll: "
## named integer(0)
## [1] "numeric data w/ NaNs in glbObsAll: "
## named integer(0)
## [1] "string data missing in glbObsAll: "
##       NewsDesk    SectionName SubsectionName       Headline        Snippet 
##           2408           2899           6176              0             13 
##       Abstract        PubDate 
##             17              0
```

```
##   Popular Pplr.fctr   .n
## 1       0         N 5439
## 2      NA      <NA> 1870
## 3       1         Y 1093
```

```
## Warning: Removed 1 rows containing missing values (position_stack).
```

![](NYTBlogs3_feat_PubDate_files/figure-html/inspect.data-2.png) 

```
##       Pplr.fctr.N Pplr.fctr.Y Pplr.fctr.NA
## Test           NA          NA         1870
## Train        5439        1093           NA
##       Pplr.fctr.N Pplr.fctr.Y Pplr.fctr.NA
## Test           NA          NA            1
## Train   0.8326699   0.1673301           NA
```

```
## Loading required package: plyr
## -------------------------------------------------------------------------
## You have loaded plyr after dplyr - this is likely to cause problems.
## If you need functions from both plyr and dplyr, please load plyr first, then dplyr:
## library(plyr); library(dplyr)
## -------------------------------------------------------------------------
## 
## Attaching package: 'plyr'
## 
## The following objects are masked from 'package:dplyr':
## 
##     arrange, count, desc, failwith, id, mutate, rename, summarise,
##     summarize
## 
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](NYTBlogs3_feat_PubDate_files/figure-html/inspect.data-3.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/inspect.data-4.png) 

```
##          label step_major step_minor label_minor    bgn    end elapsed
## 2 inspect.data          2          0           0 21.011 25.368   4.357
## 3   scrub.data          2          1           1 25.368     NA      NA
```

### Step `2.1: scrub data`

```
## [1] "numeric data missing in glbObsAll: "
##   Popular Pplr.fctr 
##      1870      1870 
## [1] "numeric data w/ 0s in glbObsAll: "
## WordCount   Popular 
##       109      5439 
## [1] "numeric data w/ Infs in glbObsAll: "
## named integer(0)
## [1] "numeric data w/ NaNs in glbObsAll: "
## named integer(0)
## [1] "string data missing in glbObsAll: "
##       NewsDesk    SectionName SubsectionName       Headline        Snippet 
##           2408           2899           6176              0             13 
##       Abstract        PubDate 
##             17              0
```

```
##            label step_major step_minor label_minor    bgn    end elapsed
## 3     scrub.data          2          1           1 25.368 26.528    1.16
## 4 transform.data          2          2           2 26.529     NA      NA
```

### Step `2.2: transform data`

```
## [1] "Creating new feature: NDSS.my..."
## [1] "Creating new feature: WordCount.log1p..."
## [1] "Creating new feature: WordCount.root2..."
## [1] "Creating new feature: WordCount.nexp..."
```

```
##              label step_major step_minor label_minor    bgn    end elapsed
## 4   transform.data          2          2           2 26.529 26.626   0.097
## 5 extract.features          3          0           0 26.626     NA      NA
```

## Step `3.0: extract features`

```
##                  label step_major step_minor label_minor    bgn end
## 1 extract.features_bgn          1          0           0 26.682  NA
##   elapsed
## 1      NA
```

```
##                                 label step_major step_minor label_minor
## 1                extract.features_bgn          1          0           0
## 2 extract.features_factorize.str.vars          2          0           0
##      bgn    end elapsed
## 1 26.682 26.701   0.019
## 2 26.701     NA      NA
```

```
##         NewsDesk      SectionName   SubsectionName         Headline 
##       "NewsDesk"    "SectionName" "SubsectionName"       "Headline" 
##          Snippet         Abstract          PubDate             .src 
##        "Snippet"       "Abstract"        "PubDate"           ".src" 
##          NDSS.my 
##        "NDSS.my"
```

```
## Warning: Creating factors of string variable: NDSS.my: # of unique values:
## 21
```

```
##                                   label step_major step_minor label_minor
## 2   extract.features_factorize.str.vars          2          0           0
## 3 extract.features_xtract.DateTime.vars          3          0           0
##      bgn   end elapsed
## 2 26.701 26.72   0.019
## 3 26.720    NA      NA
## [1] "Extracting features from DateTime(s): PubDate"
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](NYTBlogs3_feat_PubDate_files/figure-html/extract.features-1.png) 

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
## Loading required package: XML
```

![](NYTBlogs3_feat_PubDate_files/figure-html/extract.features-2.png) 

```
## [1] "**********"
## [1] "Consider adding state & city holidays for glbFeatsDateTime: PubDate"
## [1] "**********"
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](NYTBlogs3_feat_PubDate_files/figure-html/extract.features-3.png) 

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](NYTBlogs3_feat_PubDate_files/figure-html/extract.features-4.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/extract.features-5.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/extract.features-6.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/extract.features-7.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/extract.features-8.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/extract.features-9.png) 

```
## Loading required package: zoo
## 
## Attaching package: 'zoo'
## 
## The following objects are masked from 'package:base':
## 
##     as.Date, as.Date.numeric
```

![](NYTBlogs3_feat_PubDate_files/figure-html/extract.features-10.png) 

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](NYTBlogs3_feat_PubDate_files/figure-html/extract.features-11.png) 

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](NYTBlogs3_feat_PubDate_files/figure-html/extract.features-12.png) 

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](NYTBlogs3_feat_PubDate_files/figure-html/extract.features-13.png) 

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](NYTBlogs3_feat_PubDate_files/figure-html/extract.features-14.png) 

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](NYTBlogs3_feat_PubDate_files/figure-html/extract.features-15.png) 

```
##                                   label step_major step_minor label_minor
## 3 extract.features_xtract.DateTime.vars          3          0           0
## 4                  extract.features_end          4          0           0
##      bgn    end elapsed
## 3 26.720 48.247  21.527
## 4 48.247     NA      NA
```

```
##                                   label step_major step_minor label_minor
## 3 extract.features_xtract.DateTime.vars          3          0           0
## 1                  extract.features_bgn          1          0           0
## 2   extract.features_factorize.str.vars          2          0           0
##      bgn    end elapsed duration
## 3 26.720 48.247  21.527   21.527
## 1 26.682 26.701   0.019    0.019
## 2 26.701 26.720   0.019    0.019
## [1] "Total Elapsed Time: 48.247 secs"
```

![](NYTBlogs3_feat_PubDate_files/figure-html/extract.features-16.png) 

```
## time	trans	 "bgn " "fit.data.training.all " "predict.data.new " "end " 
## 0.0000 	multiple enabled transitions:  data.training.all data.new model.selected 	firing:  data.training.all 
## 1.0000 	 1 	 2 1 0 0 
## 1.0000 	multiple enabled transitions:  data.training.all data.new model.selected model.final data.training.all.prediction 	firing:  data.new 
## 2.0000 	 2 	 1 1 1 0
```

![](NYTBlogs3_feat_PubDate_files/figure-html/extract.features-17.png) 

```
##                 label step_major step_minor label_minor    bgn    end
## 5    extract.features          3          0           0 26.626 49.557
## 6 manage.missing.data          3          1           1 49.558     NA
##   elapsed
## 5  22.932
## 6      NA
```

### Step `3.1: manage missing data`

```
## [1] "numeric data missing in : "
##   Popular Pplr.fctr 
##      1870      1870 
## [1] "numeric data w/ 0s in : "
##            WordCount              Popular      WordCount.log1p 
##                  109                 5439                  109 
##      WordCount.root2       WordCount.nexp   PubDate.wkday.fctr 
##                  109                 2044                  378 
##        PubDate.wkend        PubDate.hlday  PubDate.day.minutes 
##                 7787                 8160                    5 
##  PubDate.last2.log1p  PubDate.last4.log1p  PubDate.last8.log1p 
##                    2                    4                    8 
## PubDate.last16.log1p PubDate.last32.log1p 
##                   16                   32 
## [1] "numeric data w/ Infs in : "
## named integer(0)
## [1] "numeric data w/ NaNs in : "
## named integer(0)
## [1] "string data missing in : "
##       NewsDesk    SectionName SubsectionName       Headline        Snippet 
##           2408           2899           6176              0             13 
##       Abstract        PubDate        NDSS.my 
##             17              0              0
```

```
## [1] "numeric data missing in : "
##   Popular Pplr.fctr 
##      1870      1870 
## [1] "numeric data w/ 0s in : "
##            WordCount              Popular      WordCount.log1p 
##                  109                 5439                  109 
##      WordCount.root2       WordCount.nexp   PubDate.wkday.fctr 
##                  109                 2044                  378 
##        PubDate.wkend        PubDate.hlday  PubDate.day.minutes 
##                 7787                 8160                    5 
##  PubDate.last2.log1p  PubDate.last4.log1p  PubDate.last8.log1p 
##                    2                    4                    8 
## PubDate.last16.log1p PubDate.last32.log1p 
##                   16                   32 
## [1] "numeric data w/ Infs in : "
## named integer(0)
## [1] "numeric data w/ NaNs in : "
## named integer(0)
## [1] "string data missing in : "
##       NewsDesk    SectionName SubsectionName       Headline        Snippet 
##           2408           2899           6176              0             13 
##       Abstract        PubDate        NDSS.my 
##             17              0              0
```

```
##                 label step_major step_minor label_minor    bgn    end
## 6 manage.missing.data          3          1           1 49.558 56.087
## 7        cluster.data          3          2           2 56.087     NA
##   elapsed
## 6   6.529
## 7      NA
```

## Step `3.2: cluster data`

```r
mycompute_entropy_df <- function(obs_df, entropy_var, by_var=NULL) {   
    require(lazyeval)
    require(dplyr)
    require(tidyr)

    if (is.null(by_var)) {
        by_var <- ".default"
        obs_df$.default <- as.factor(".default") 
    }
    
    if (!any(grepl(".clusterid", names(obs_df), fixed=TRUE)))
        obs_df$.clusterid <- 1
        
    cluster_df <- obs_df %>%
            count_(c(by_var, ".clusterid", entropy_var)) %>%
            dplyr::filter(n > 0) %>%
            dplyr::filter_(interp(~(!is.na(var)), var=as.name(entropy_var))) %>%
            unite_(paste0(by_var, ".clusterid"),
                   c(interp(by_var), ".clusterid")) %>%
            spread_(interp(entropy_var), "n", fill=0) 

#     head(cluster_df)
#     sum(cluster_df$n)
    tmp.entropy <- sapply(1:nrow(cluster_df),
            function(row) entropy(as.numeric(cluster_df[row, -1]), method = "ML"))
    tmp.knt <- sapply(1:nrow(cluster_df),
                    function(row) sum(as.numeric(cluster_df[row, -1])))
    cluster_df$.entropy <- tmp.entropy; cluster_df$.knt <- tmp.knt
    #print(cluster_df)
    return(cluster_df)
}
    
if (glb_cluster) {
    require(proxy)
    #require(hash)
    require(dynamicTreeCut)
    require(entropy)
    require(tidyr)
    require(ggdendro)

    mywgtdcosine_dist <- function(x, y=NULL, weights=NULL) {
        if (!inherits(x, "matrix"))
            x <- as.matrix(x)
    
        if (is.null(weights))
            weights <- rep(1, ncol(x))
    
        wgtsx <- matrix(rep(weights / sum(weights), nrow(x)), nrow = nrow(x),
                        byrow = TRUE)
        wgtdx <- x * wgtsx
    
        wgtdxsqsum <- as.matrix(rowSums((x ^ 2) * wgtsx), byrow=FALSE)
        denom <- sqrt(wgtdxsqsum %*% t(wgtdxsqsum))
    
        ret_mtrx <- 1 - ((sum(weights) ^ 1) * (wgtdx %*% t(wgtdx)) / denom)
        ret_mtrx[is.nan(ret_mtrx)] <- 1
        diag(ret_mtrx) <- 0
        return(ret_mtrx)
    }
    #pr_DB$delete_entry("mywgtdcosine"); 
    # Need to do this only once across runs ?
    if (!pr_DB$entry_exists("mywgtdcosine")) {
        pr_DB$set_entry(FUN = mywgtdcosine_dist, names = c("mywgtdcosine"))
        pr_DB$modify_entry(names="mywgtdcosine", type="metric", loop=FALSE)
    }
    #pr_DB$get_entry("mywgtdcosine")

#     glb_hash <- hash(key=unique(glbObsAll$myCategory), 
#                      values=1:length(unique(glbObsAll$myCategory)))
#     glb_hash_lst <- hash(key=unique(glbObsAll$myCategory), 
#                      values=1:length(unique(glbObsAll$myCategory)))
#stop(here"); glb2Sav(); glbObsAll <- savObsAll
    cluster_vars <- grep(paste0("[", 
                        toupper(paste0(substr(glbFeatsText, 1, 1), collapse = "")),
                                      "]\\.[PT]\\."), 
                               names(glbObsAll), value = TRUE)
    # Assign correlations with rsp_var as weights for cosine distance
    print("Clustering features: ")
    cluster_vars_df <- data.frame(abs.cor.y = abs(cor(
                        glbObsAll[glbObsAll$.src == "Train", cluster_vars],
            as.numeric(glbObsAll[glbObsAll$.src == "Train", glb_rsp_var]),
                                    use = "pairwise.complete.obs")))
    print(tail(cluster_vars_df <- orderBy(~ abs.cor.y, 
                                    subset(cluster_vars_df, !is.na(abs.cor.y))), 5))
    print(sprintf("    .rnorm cor: %0.4f",
        cor(glbObsAll[glbObsAll$.src == "Train", ".rnorm"], 
            as.numeric(glbObsAll[glbObsAll$.src == "Train", glb_rsp_var]),
            use = "pairwise.complete.obs")))
    
    print(sprintf("glbObsAll Entropy: %0.4f", 
        allobs_ent <- entropy(table(glbObsAll[, glb_cluster_entropy_var]),
                              method="ML")))
    
    print(category_df <- mycompute_entropy_df(obs_df=glbObsAll,
                                             entropy_var=glb_cluster_entropy_var,
                                             by_var=glbFeatsCategory))
    print(sprintf("glbObsAll$%s Entropy: %0.4f (%0.4f pct)",
                    glbFeatsCategory,
            category_ent <- weighted.mean(category_df$.entropy, category_df$.knt),
                    100 * category_ent / allobs_ent))

    glbObsAll$.clusterid <- 1    
    #print(max(table(glbObsAll$myCategory.fctr) / 20))
        
#stop(here"); glb2Sav()    
    grp_ids <- sort(unique(glbObsAll[, glbFeatsCategory]))
    glb_cluster_size_df_lst <- list()
    png(paste0(glb_out_pfx, "FeatsTxtClusters.png"), 
        width = 480 * 2, height = 480 * length(grp_ids))
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow = length(grp_ids), ncol = 2)))
    pltIx <- 1
    for (grp in grp_ids) {
# if (grep(grp, levels(grp_ids)) <= 6) next                
# if (grep(grp, levels(grp_ids)) > 9) next        
# if (grep(grp, levels(grp_ids)) != 10) next        
        print(sprintf("Category: %s", grp))
        ctgry_allobs_df <- glbObsAll[glbObsAll[, glbFeatsCategory] == grp, ]
        if (!inherits(ctgry_allobs_df[, glb_cluster_entropy_var], "factor"))
            ctgry_allobs_df[, glb_cluster_entropy_var] <- 
                as.factor(ctgry_allobs_df[, glb_cluster_entropy_var])
        
        #dstns_dist <- proxy::dist(ctgry_allobs_df[, cluster_vars], method = "cosine")
        dstns_dist <- proxy::dist(ctgry_allobs_df[, row.names(cluster_vars_df)], 
                                  method = "mywgtdcosine",
                                  weights = cluster_vars_df$abs.cor.y)
        # Custom distance functions return a crossdist object
        #dstns_mtrx <- as.matrix(dstns_dist)
        dstns_mtrx <- matrix(as.vector(dstns_dist), nrow=attr(dstns_dist, "dim")[1],
                             dimnames=attr(dstns_dist, "dimnames"))
        dstns_dist <- as.dist(dstns_mtrx)

        print(sprintf("max distance(%0.4f) pair:", max(dstns_mtrx)))
#         print(dim(dstns_mtrx))        
#         print(sprintf("which.max: %d", which.max(dstns_mtrx)))
        row_ix <- ceiling(which.max(dstns_mtrx) / ncol(dstns_mtrx))
        col_ix <- which.max(dstns_mtrx[row_ix, ])
#         print(sprintf("row_ix: %d", row_ix)); print(sprintf("col_ix: %d", col_ix));
#         print(dim(ctgry_allobs_df))
        print(ctgry_allobs_df[c(row_ix, col_ix), 
            c(glb_id_var, glb_cluster_entropy_var, glbFeatsCategory, glbFeatsText, cluster_vars)])
    
        min_dstns_mtrx <- dstns_mtrx
        diag(min_dstns_mtrx) <- 1
        # Float representations issue -2.22e-16 vs. 0.0000
        print(sprintf("min distance(%0.4f) pair:", min(min_dstns_mtrx)))
        row_ix <- ceiling(which.min(min_dstns_mtrx) / ncol(min_dstns_mtrx))
        col_ix <- which.min(min_dstns_mtrx[row_ix, ])
        print(ctgry_allobs_df[c(row_ix, col_ix), 
            c(glb_id_var, glb_cluster_entropy_var, glbFeatsCategory, glbFeatsText,
              cluster_vars)])
    
        set.seed(glb_cluster.seed)
        clusters <- hclust(dstns_dist, method = "ward.D2")
        # Workaround to avoid "Error in cutree(dendro, h = heightcutoff) : the 'height' component of 'tree' is not sorted (increasingly)"
        if (with(clusters,all.equal(height,sort(height))))
            clusters$height <- round(clusters$height,6)
        
        clusters$labels <- ctgry_allobs_df[, glb_id_var]
        clustersDD <- dendro_data(clusters)
        clustersDD$labels[, glb_rsp_var] <- sapply(clustersDD$labels$label, function(id)
            ctgry_allobs_df[id == ctgry_allobs_df[, glb_id_var], glb_rsp_var])
        print(ggdendrogram(clustersDD, rotate = TRUE, size = 2) + 
                geom_point(data = clustersDD$labels, 
            aes_string(x = "x", color = glb_rsp_var), y = min(clustersDD$segments$y)) + 
                coord_flip(ylim = c(min(clustersDD$segments$y),
                                         max(clustersDD$segments$y))) + 
                ggtitle(grp),
            vp = viewport(layout.pos.row = pltIx, layout.pos.col = 1))  
        
#         clusters$labels <- ctgry_allobs_df[, glb_id_var]
#         clustersDD <- dendro_data(clusters)
#         clustersDD$labels$color <- sapply(clustersDD$labels$label, function(id)
#             ctgry_allobs_df[id == ctgry_allobs_df[, glb_id_var], glb_rsp_var])
#         print(ggdendrogram(clustersDD, rotate = TRUE, size = 2) + 
#                 geom_point(data = clustersDD$labels, 
#                 aes_string(x = "x", color = "color"), y = min(clustersDD$segments$y)) + 
#                 coord_flip(ylim = c(min(clustersDD$segments$y),
#                                          max(clustersDD$segments$y))))
#         print(ggdendrogram(clustersDD, rotate = TRUE, size = 2) + 
#                 geom_point(data = clustersDD$labels, 
#                           aes_string(x = "x", y = "y", color = "color")))
#         myplclust(clusters, lab=ctgry_allobs_df[, glb_id_var], 
#                   lab.col=unclass(ctgry_allobs_df[, glb_cluster_entropy_var]))

        opt_minclustersize_df <- data.frame(minclustersize = nrow(ctgry_allobs_df), 
            entropy = entropy(table(ctgry_allobs_df[, glb_cluster_entropy_var]),
                              method = "ML"))
        for (minclustersize in 
             as.integer(seq(nrow(ctgry_allobs_df) / 2, nrow(ctgry_allobs_df) / 10, 
                            length = 5))) {
            clusterGroups <- cutreeDynamic(clusters, minClusterSize = minclustersize,
                                           method = "tree", deepSplit = 0)
            # Unassigned groups are labeled 0; the largest group has label 1
            clusterGroups[clusterGroups == 0] <- 1
            ctgry_allobs_df$.clusterid <- clusterGroups
            ctgry_clstrs_df <- mycompute_entropy_df(ctgry_allobs_df,
                                                    glb_cluster_entropy_var)
            opt_minclustersize_df <- rbind(opt_minclustersize_df, 
                                           data.frame(minclustersize = minclustersize,
                entropy = weighted.mean(ctgry_clstrs_df$.entropy, ctgry_clstrs_df$.knt)))
        }
        opt_minclustersize <-
            opt_minclustersize_df$minclustersize[which.min(opt_minclustersize_df$entropy)]
        opt_minclustersize_df$.color <- 
            ifelse(opt_minclustersize_df$minclustersize == opt_minclustersize,
                   "red", "blue")
        print(ggplot(data = opt_minclustersize_df, 
                     mapping = aes(x = minclustersize, y = entropy)) + 
                geom_point(aes(color = .color)) + scale_color_identity() + 
                guides(color = "none") + geom_line(),
            vp = viewport(layout.pos.row = pltIx, layout.pos.col = 2))
        glb_cluster_size_df_lst[[grp]] <- opt_minclustersize_df
        
        # select minclustersize that minimizes entropy
        clusterGroups <- cutreeDynamic(clusters, minClusterSize = opt_minclustersize,
                                       method = "tree",
                                       deepSplit = 0)
        # Unassigned groups are labeled 0; the largest group has label 1
        table(clusterGroups, ctgry_allobs_df[, glb_cluster_entropy_var], 
              useNA = "ifany")   
        clusterGroups[clusterGroups == 0] <- 1
        table(clusterGroups, ctgry_allobs_df[, glb_cluster_entropy_var], useNA = "ifany") 
        glbObsAll[glbObsAll[, glbFeatsCategory] == grp,]$.clusterid <-
            clusterGroups
        
        pltIx <- pltIx + 1
    }
    dev.off()
    #all.equal(savObsAll_clusterid, glbObsAll$.clusterid)
    
    print(cluster_df <- mycompute_entropy_df(obs_df=glbObsAll,
                                             entropy_var=glb_cluster_entropy_var,
                                             by_var=glbFeatsCategory))
    print(sprintf("glbObsAll$%s$.clusterid Entropy: %0.4f (%0.4f pct)",
                    glbFeatsCategory,
                cluster_ent <- weighted.mean(cluster_df$.entropy, cluster_df$.knt),
                    100 * cluster_ent / category_ent))
    
    glbObsAll$.clusterid.fctr <- as.factor(glbObsAll$.clusterid)
    # .clusterid.fctr is created automatically (probably ?) later
    glbFeatsExclude <- c(glbFeatsExclude, ".clusterid")
    if (!is.null(glbFeatsCategory))
#         glbFeatsInteractionOnly[ifelse(grepl("\\.fctr", glbFeatsCategory),
#                                             glbFeatsCategory, 
#                                             paste0(glbFeatsCategory, ".fctr"))] <-
#             c(".clusterid.fctr")
        glbFeatsInteractionOnly[[".clusterid.fctr"]] <-
            ifelse(grepl("\\.fctr", glbFeatsCategory), glbFeatsCategory, 
                                                        paste0(glbFeatsCategory, ".fctr"))
            
    if (glbFeatsTextClusterVarsExclude)
        glbFeatsExclude <- c(glbFeatsExclude, cluster_vars)
}

# Last call for data modifications 
#stop(here") # savObsAll <- glbObsAll
# glbObsAll[(glbObsAll$PropR == 0.75) & (glbObsAll$State == "Hawaii"), "PropR.fctr"] <- "N"

# Re-partition
glbObsTrn <- subset(glbObsAll, .src == "Train")
glbObsNew <- subset(glbObsAll, .src == "Test")

glb_chunks_df <- myadd_chunk(glb_chunks_df, "partition.data.training", major.inc = TRUE)
```

```
##                     label step_major step_minor label_minor    bgn    end
## 7            cluster.data          3          2           2 56.087 56.151
## 8 partition.data.training          4          0           0 56.152     NA
##   elapsed
## 7   0.064
## 8      NA
```

## Step `4.0: partition data training`

```
## [1] "Prediction Hints by Catgeory:"
##    NDSS.my.fctr Popular.0 Popular.1 .n.tst .strata.0 .strata.1
## 5   #U.S.#Edctn       325        NA     89        82        17
## 10       Cltr##         1        NA     70         1        13
## 12   Frgn#Wrld#       172        NA     47        44         9
## 21       myOthr        38        NA      5         5         1
```

```
## Loading required package: sampling
## 
## Attaching package: 'sampling'
## 
## The following objects are masked from 'package:survival':
## 
##     cluster, strata
## 
## The following object is masked from 'package:caret':
## 
##     cluster
```

```
##     Popular.0 Popular.1 Popular.NA
##            NA        NA       1870
## Fit      3941       863         NA
## OOB      1498       230         NA
##     Popular.0 Popular.1 Popular.NA
##            NA        NA          1
## Fit 0.8203580 0.1796420         NA
## OOB 0.8668981 0.1331019         NA
##               NDSS.my.fctr .n.Fit .n.OOB .n.Tst .freqRatio.Fit
## 1                       ##    913    371    342    0.190049958
## 6       Bsnss#BsnssDy#Dlbk    629    323    304    0.130932556
## 11              Cltr#Arts#    490    185    174    0.101998335
## 15              OpEd#Opnn#    437     89    164    0.090965862
## 9           Bsnss#Tchnlgy#    213    126    114    0.044338052
## 19                 TStyl##    623    101    105    0.129683597
## 5              #U.S.#Edctn    243     82     89    0.050582848
## 10                  Cltr##     NA      1     70             NA
## 14           Mtr#N.Y./Rgn#    128     70     67    0.026644463
## 18             Styls#U.S.#    127     50     61    0.026436303
## 16              Scnc#Hlth#    148     48     57    0.030807660
## 13        Frgn#Wrld#AsPcfc    150     53     56    0.031223980
## 2                  #Mltmd#     92     49     52    0.019150708
## 12              Frgn#Wrld#    128     44     47    0.026644463
## 8      Bsnss#Crsswrds/Gms#    105     18     42    0.021856786
## 7  Bsnss#BsnssDy#SmllBsnss    100     40     41    0.020815987
## 20              Trvl#Trvl#     83     34     35    0.017277269
## 3            #Opnn#RmFrDbt     42     20     20    0.008742714
## 17             Styls##Fshn    104     15     15    0.021648626
## 4         #Opnn#ThPblcEdtr     16      4     10    0.003330558
## 21                  myOthr     33      5      5    0.006869276
##    .freqRatio.OOB .freqRatio.Tst
## 1    0.2146990741    0.182887701
## 6    0.1869212963    0.162566845
## 11   0.1070601852    0.093048128
## 15   0.0515046296    0.087700535
## 9    0.0729166667    0.060962567
## 19   0.0584490741    0.056149733
## 5    0.0474537037    0.047593583
## 10   0.0005787037    0.037433155
## 14   0.0405092593    0.035828877
## 18   0.0289351852    0.032620321
## 16   0.0277777778    0.030481283
## 13   0.0306712963    0.029946524
## 2    0.0283564815    0.027807487
## 12   0.0254629630    0.025133690
## 8    0.0104166667    0.022459893
## 7    0.0231481481    0.021925134
## 20   0.0196759259    0.018716578
## 3    0.0115740741    0.010695187
## 17   0.0086805556    0.008021390
## 4    0.0023148148    0.005347594
## 21   0.0028935185    0.002673797
```

```
## [1] "glbObsAll: "
```

```
## [1] 8402   42
```

```
## [1] "glbObsTrn: "
```

```
## [1] 6532   42
```

```
## [1] "glbObsFit: "
```

```
## [1] 4804   41
```

```
## [1] "glbObsOOB: "
```

```
## [1] 1728   41
```

```
## [1] "glbObsNew: "
```

```
## [1] 1870   41
```

```
## Warning in rm(split): object 'split' not found
```

```
##                     label step_major step_minor label_minor    bgn    end
## 8 partition.data.training          4          0           0 56.152 57.658
## 9         select.features          5          0           0 57.659     NA
##   elapsed
## 8   1.507
## 9      NA
```

## Step `5.0: select features`

```
## Warning in cor(data.matrix(entity_df[, sel_feats]), y =
## as.numeric(entity_df[, : the standard deviation is zero
```

```
##                                                    id        cor.y
## Popular                                       Popular  1.000000000
## WordCount.root2                       WordCount.root2  0.292120679
## WordCount                                   WordCount  0.257526549
## WordCount.log1p                       WordCount.log1p  0.254319628
## NDSS.my.fctr                             NDSS.my.fctr  0.165445970
## PubDate.day.minutes               PubDate.day.minutes  0.156753478
## PubDate.day.minutes.poly.1 PubDate.day.minutes.poly.1  0.156753478
## PubDate.hour.fctr                   PubDate.hour.fctr  0.135436805
## PubDate.wkend                           PubDate.wkend  0.104707290
## PubDate.day.minutes.poly.4 PubDate.day.minutes.poly.4  0.073941394
## PubDate.day.minutes.poly.2 PubDate.day.minutes.poly.2  0.070977720
## PubDate.last4.log1p               PubDate.last4.log1p  0.066473282
## PubDate.last2.log1p               PubDate.last2.log1p  0.063068716
## PubDate.day.minutes.poly.5 PubDate.day.minutes.poly.5 -0.055929231
## PubDate.last8.log1p               PubDate.last8.log1p  0.054458821
## WordCount.nexp                         WordCount.nexp -0.053208396
## PubDate.last16.log1p             PubDate.last16.log1p  0.040735543
## PubDate.wkday.fctr                 PubDate.wkday.fctr -0.039801288
## PubDate.minute.fctr               PubDate.minute.fctr -0.034073846
## PubDate.day.minutes.poly.3 PubDate.day.minutes.poly.3  0.027983551
## PubDate.month.fctr                 PubDate.month.fctr  0.019148739
## PubDate.POSIX                           PubDate.POSIX  0.015683258
## PubDate.hlday                           PubDate.hlday  0.014690122
## PubDate.juliandate                 PubDate.juliandate  0.014361075
## PubDate.zoo                               PubDate.zoo  0.013260902
## PubDate.second.fctr               PubDate.second.fctr -0.011879458
## UniqueID                                     UniqueID  0.011824920
## PubDate.date.fctr                   PubDate.date.fctr -0.011647558
## .rnorm                                         .rnorm  0.008212201
## PubDate.last32.log1p             PubDate.last32.log1p  0.003558081
## PubDate.year.fctr                   PubDate.year.fctr           NA
##                            exclude.as.feat   cor.y.abs
## Popular                                  1 1.000000000
## WordCount.root2                          0 0.292120679
## WordCount                                1 0.257526549
## WordCount.log1p                          0 0.254319628
## NDSS.my.fctr                             0 0.165445970
## PubDate.day.minutes                      1 0.156753478
## PubDate.day.minutes.poly.1               0 0.156753478
## PubDate.hour.fctr                        0 0.135436805
## PubDate.wkend                            0 0.104707290
## PubDate.day.minutes.poly.4               0 0.073941394
## PubDate.day.minutes.poly.2               0 0.070977720
## PubDate.last4.log1p                      0 0.066473282
## PubDate.last2.log1p                      0 0.063068716
## PubDate.day.minutes.poly.5               0 0.055929231
## PubDate.last8.log1p                      0 0.054458821
## WordCount.nexp                           0 0.053208396
## PubDate.last16.log1p                     0 0.040735543
## PubDate.wkday.fctr                       0 0.039801288
## PubDate.minute.fctr                      0 0.034073846
## PubDate.day.minutes.poly.3               0 0.027983551
## PubDate.month.fctr                       0 0.019148739
## PubDate.POSIX                            1 0.015683258
## PubDate.hlday                            0 0.014690122
## PubDate.juliandate                       0 0.014361075
## PubDate.zoo                              1 0.013260902
## PubDate.second.fctr                      0 0.011879458
## UniqueID                                 1 0.011824920
## PubDate.date.fctr                        0 0.011647558
## .rnorm                                   0 0.008212201
## PubDate.last32.log1p                     0 0.003558081
## PubDate.year.fctr                        0          NA
```

```
## [1] "cor(PubDate.juliandate, PubDate.month.fctr)=0.9393"
## [1] "cor(Pplr.fctr, PubDate.juliandate)=0.0144"
## [1] "cor(Pplr.fctr, PubDate.month.fctr)=0.0191"
```

```
## Warning in myfind_cor_features(feats_df = glb_feats_df, obs_df =
## glbObsTrn, : Identified PubDate.juliandate as highly correlated with
## PubDate.month.fctr
```

```
## [1] "cor(PubDate.day.minutes.poly.1, PubDate.hour.fctr)=0.9026"
## [1] "cor(Pplr.fctr, PubDate.day.minutes.poly.1)=0.1568"
## [1] "cor(Pplr.fctr, PubDate.hour.fctr)=0.1354"
```

```
## Warning in myfind_cor_features(feats_df = glb_feats_df, obs_df =
## glbObsTrn, : Identified PubDate.hour.fctr as highly correlated with
## PubDate.day.minutes.poly.1
```

```
## [1] "cor(WordCount.log1p, WordCount.root2)=0.8906"
## [1] "cor(Pplr.fctr, WordCount.log1p)=0.2543"
## [1] "cor(Pplr.fctr, WordCount.root2)=0.2921"
```

```
## Warning in myfind_cor_features(feats_df = glb_feats_df, obs_df =
## glbObsTrn, : Identified WordCount.log1p as highly correlated with
## WordCount.root2
```

```
## [1] "cor(PubDate.last4.log1p, PubDate.last8.log1p)=0.8253"
## [1] "cor(Pplr.fctr, PubDate.last4.log1p)=0.0665"
## [1] "cor(Pplr.fctr, PubDate.last8.log1p)=0.0545"
```

```
## Warning in myfind_cor_features(feats_df = glb_feats_df, obs_df =
## glbObsTrn, : Identified PubDate.last8.log1p as highly correlated with
## PubDate.last4.log1p
```

```
## [1] "cor(PubDate.last2.log1p, PubDate.last4.log1p)=0.7598"
## [1] "cor(Pplr.fctr, PubDate.last2.log1p)=0.0631"
## [1] "cor(Pplr.fctr, PubDate.last4.log1p)=0.0665"
```

```
## Warning in myfind_cor_features(feats_df = glb_feats_df, obs_df =
## glbObsTrn, : Identified PubDate.last2.log1p as highly correlated with
## PubDate.last4.log1p
```

```
##                                                    id        cor.y
## Popular                                       Popular  1.000000000
## WordCount.root2                       WordCount.root2  0.292120679
## WordCount                                   WordCount  0.257526549
## WordCount.log1p                       WordCount.log1p  0.254319628
## NDSS.my.fctr                             NDSS.my.fctr  0.165445970
## PubDate.day.minutes               PubDate.day.minutes  0.156753478
## PubDate.day.minutes.poly.1 PubDate.day.minutes.poly.1  0.156753478
## PubDate.hour.fctr                   PubDate.hour.fctr  0.135436805
## PubDate.wkend                           PubDate.wkend  0.104707290
## PubDate.day.minutes.poly.4 PubDate.day.minutes.poly.4  0.073941394
## PubDate.day.minutes.poly.2 PubDate.day.minutes.poly.2  0.070977720
## PubDate.last4.log1p               PubDate.last4.log1p  0.066473282
## PubDate.last2.log1p               PubDate.last2.log1p  0.063068716
## PubDate.last8.log1p               PubDate.last8.log1p  0.054458821
## PubDate.last16.log1p             PubDate.last16.log1p  0.040735543
## PubDate.day.minutes.poly.3 PubDate.day.minutes.poly.3  0.027983551
## PubDate.month.fctr                 PubDate.month.fctr  0.019148739
## PubDate.POSIX                           PubDate.POSIX  0.015683258
## PubDate.hlday                           PubDate.hlday  0.014690122
## PubDate.juliandate                 PubDate.juliandate  0.014361075
## PubDate.zoo                               PubDate.zoo  0.013260902
## UniqueID                                     UniqueID  0.011824920
## .rnorm                                         .rnorm  0.008212201
## PubDate.last32.log1p             PubDate.last32.log1p  0.003558081
## PubDate.date.fctr                   PubDate.date.fctr -0.011647558
## PubDate.second.fctr               PubDate.second.fctr -0.011879458
## PubDate.minute.fctr               PubDate.minute.fctr -0.034073846
## PubDate.wkday.fctr                 PubDate.wkday.fctr -0.039801288
## WordCount.nexp                         WordCount.nexp -0.053208396
## PubDate.day.minutes.poly.5 PubDate.day.minutes.poly.5 -0.055929231
## PubDate.year.fctr                   PubDate.year.fctr           NA
##                            exclude.as.feat   cor.y.abs
## Popular                                  1 1.000000000
## WordCount.root2                          0 0.292120679
## WordCount                                1 0.257526549
## WordCount.log1p                          0 0.254319628
## NDSS.my.fctr                             0 0.165445970
## PubDate.day.minutes                      1 0.156753478
## PubDate.day.minutes.poly.1               0 0.156753478
## PubDate.hour.fctr                        0 0.135436805
## PubDate.wkend                            0 0.104707290
## PubDate.day.minutes.poly.4               0 0.073941394
## PubDate.day.minutes.poly.2               0 0.070977720
## PubDate.last4.log1p                      0 0.066473282
## PubDate.last2.log1p                      0 0.063068716
## PubDate.last8.log1p                      0 0.054458821
## PubDate.last16.log1p                     0 0.040735543
## PubDate.day.minutes.poly.3               0 0.027983551
## PubDate.month.fctr                       0 0.019148739
## PubDate.POSIX                            1 0.015683258
## PubDate.hlday                            0 0.014690122
## PubDate.juliandate                       0 0.014361075
## PubDate.zoo                              1 0.013260902
## UniqueID                                 1 0.011824920
## .rnorm                                   0 0.008212201
## PubDate.last32.log1p                     0 0.003558081
## PubDate.date.fctr                        0 0.011647558
## PubDate.second.fctr                      0 0.011879458
## PubDate.minute.fctr                      0 0.034073846
## PubDate.wkday.fctr                       0 0.039801288
## WordCount.nexp                           0 0.053208396
## PubDate.day.minutes.poly.5               0 0.055929231
## PubDate.year.fctr                        0          NA
##                                            cor.high.X freqRatio
## Popular                                          <NA>  4.976212
## WordCount.root2                                  <NA>  2.315789
## WordCount                                        <NA>  2.315789
## WordCount.log1p                       WordCount.root2  2.315789
## NDSS.my.fctr                                     <NA>  1.348739
## PubDate.day.minutes                              <NA>  1.225490
## PubDate.day.minutes.poly.1                       <NA>  1.225490
## PubDate.hour.fctr          PubDate.day.minutes.poly.1  1.835040
## PubDate.wkend                                    <NA> 12.011952
## PubDate.day.minutes.poly.4                       <NA>  1.225490
## PubDate.day.minutes.poly.2                       <NA>  1.225490
## PubDate.last4.log1p                              <NA>  1.125000
## PubDate.last2.log1p               PubDate.last4.log1p  1.375000
## PubDate.last8.log1p               PubDate.last4.log1p  1.142857
## PubDate.last16.log1p                             <NA>  3.200000
## PubDate.day.minutes.poly.3                       <NA>  1.225490
## PubDate.month.fctr                               <NA>  1.017514
## PubDate.POSIX                                    <NA>  1.000000
## PubDate.hlday                                    <NA> 28.160714
## PubDate.juliandate                 PubDate.month.fctr  1.032520
## PubDate.zoo                                      <NA>  1.000000
## UniqueID                                         <NA>  1.000000
## .rnorm                                           <NA>  1.000000
## PubDate.last32.log1p                             <NA>  8.000000
## PubDate.date.fctr                                <NA>  1.021394
## PubDate.second.fctr                              <NA>  1.018204
## PubDate.minute.fctr                              <NA>  1.483365
## PubDate.wkday.fctr                               <NA>  1.003268
## WordCount.nexp                                   <NA> 17.761364
## PubDate.day.minutes.poly.5                       <NA>  1.225490
## PubDate.year.fctr                                <NA>  0.000000
##                            percentUnique zeroVar   nzv is.cor.y.abs.low
## Popular                       0.03061849   FALSE FALSE            FALSE
## WordCount.root2              24.15799143   FALSE FALSE            FALSE
## WordCount                    24.15799143   FALSE FALSE            FALSE
## WordCount.log1p              24.15799143   FALSE FALSE            FALSE
## NDSS.my.fctr                  0.32149418   FALSE FALSE            FALSE
## PubDate.day.minutes          18.08022045   FALSE FALSE            FALSE
## PubDate.day.minutes.poly.1   18.08022045   FALSE FALSE            FALSE
## PubDate.hour.fctr             0.04592774   FALSE FALSE            FALSE
## PubDate.wkend                 0.03061849   FALSE FALSE            FALSE
## PubDate.day.minutes.poly.4   18.08022045   FALSE FALSE            FALSE
## PubDate.day.minutes.poly.2   18.08022045   FALSE FALSE            FALSE
## PubDate.last4.log1p          64.98775260   FALSE FALSE            FALSE
## PubDate.last2.log1p          51.17881200   FALSE FALSE            FALSE
## PubDate.last8.log1p          75.12247397   FALSE FALSE            FALSE
## PubDate.last16.log1p         84.44580527   FALSE FALSE            FALSE
## PubDate.day.minutes.poly.3   18.08022045   FALSE FALSE            FALSE
## PubDate.month.fctr            0.04592774   FALSE FALSE            FALSE
## PubDate.POSIX                99.86221678   FALSE FALSE            FALSE
## PubDate.hlday                 0.03061849   FALSE  TRUE            FALSE
## PubDate.juliandate            1.39314146   FALSE FALSE            FALSE
## PubDate.zoo                  99.86221678   FALSE FALSE            FALSE
## UniqueID                    100.00000000   FALSE FALSE            FALSE
## .rnorm                      100.00000000   FALSE FALSE            FALSE
## PubDate.last32.log1p         90.99816289   FALSE FALSE             TRUE
## PubDate.date.fctr             0.07654623   FALSE FALSE            FALSE
## PubDate.second.fctr           0.06123699   FALSE FALSE            FALSE
## PubDate.minute.fctr           0.06123699   FALSE FALSE            FALSE
## PubDate.wkday.fctr            0.10716473   FALSE FALSE            FALSE
## WordCount.nexp               11.32884262   FALSE FALSE            FALSE
## PubDate.day.minutes.poly.5   18.08022045   FALSE FALSE            FALSE
## PubDate.year.fctr             0.01530925    TRUE  TRUE               NA
```

```
## Warning in myplot_scatter(plt_feats_df, "percentUnique", "freqRatio",
## colorcol_name = "nzv", : converting nzv to class:factor
```

```
## Warning: Removed 9 rows containing missing values (geom_point).
```

```
## Warning: Removed 9 rows containing missing values (geom_point).
```

```
## Warning: Removed 9 rows containing missing values (geom_point).
```

![](NYTBlogs3_feat_PubDate_files/figure-html/select.features-1.png) 

```
##                                  id      cor.y exclude.as.feat  cor.y.abs
## PubDate.hlday         PubDate.hlday 0.01469012               0 0.01469012
## PubDate.year.fctr PubDate.year.fctr         NA               0         NA
##                   cor.high.X freqRatio percentUnique zeroVar  nzv
## PubDate.hlday           <NA>  28.16071    0.03061849   FALSE TRUE
## PubDate.year.fctr       <NA>   0.00000    0.01530925    TRUE TRUE
##                   is.cor.y.abs.low
## PubDate.hlday                FALSE
## PubDate.year.fctr               NA
```

![](NYTBlogs3_feat_PubDate_files/figure-html/select.features-2.png) 

```
## +(rfe) fit Fold1.Rep1 size: 55
```

```
## Warning in lda.default(x, grouping, ...): variables are collinear
```

```
## -(rfe) fit Fold1.Rep1 size: 55 
## +(rfe) imp Fold1.Rep1 
## -(rfe) imp Fold1.Rep1 
## +(rfe) fit Fold1.Rep1 size: 32 
## -(rfe) fit Fold1.Rep1 size: 32 
## +(rfe) fit Fold1.Rep1 size: 16 
## -(rfe) fit Fold1.Rep1 size: 16 
## +(rfe) fit Fold1.Rep1 size:  8 
## -(rfe) fit Fold1.Rep1 size:  8 
## +(rfe) fit Fold1.Rep1 size:  4 
## -(rfe) fit Fold1.Rep1 size:  4 
## +(rfe) fit Fold1.Rep1 size:  2 
## -(rfe) fit Fold1.Rep1 size:  2 
## +(rfe) fit Fold2.Rep1 size: 55
```

```
## Warning in lda.default(x, grouping, ...): variables are collinear
```

```
## -(rfe) fit Fold2.Rep1 size: 55 
## +(rfe) imp Fold2.Rep1 
## -(rfe) imp Fold2.Rep1 
## +(rfe) fit Fold2.Rep1 size: 32 
## -(rfe) fit Fold2.Rep1 size: 32 
## +(rfe) fit Fold2.Rep1 size: 16 
## -(rfe) fit Fold2.Rep1 size: 16 
## +(rfe) fit Fold2.Rep1 size:  8 
## -(rfe) fit Fold2.Rep1 size:  8 
## +(rfe) fit Fold2.Rep1 size:  4 
## -(rfe) fit Fold2.Rep1 size:  4 
## +(rfe) fit Fold2.Rep1 size:  2 
## -(rfe) fit Fold2.Rep1 size:  2 
## +(rfe) fit Fold3.Rep1 size: 55
```

```
## Warning in lda.default(x, grouping, ...): variables are collinear
```

```
## -(rfe) fit Fold3.Rep1 size: 55 
## +(rfe) imp Fold3.Rep1 
## -(rfe) imp Fold3.Rep1 
## +(rfe) fit Fold3.Rep1 size: 32 
## -(rfe) fit Fold3.Rep1 size: 32 
## +(rfe) fit Fold3.Rep1 size: 16 
## -(rfe) fit Fold3.Rep1 size: 16 
## +(rfe) fit Fold3.Rep1 size:  8 
## -(rfe) fit Fold3.Rep1 size:  8 
## +(rfe) fit Fold3.Rep1 size:  4 
## -(rfe) fit Fold3.Rep1 size:  4 
## +(rfe) fit Fold3.Rep1 size:  2 
## -(rfe) fit Fold3.Rep1 size:  2 
## +(rfe) fit Fold1.Rep2 size: 55
```

```
## Warning in lda.default(x, grouping, ...): variables are collinear
```

```
## -(rfe) fit Fold1.Rep2 size: 55 
## +(rfe) imp Fold1.Rep2 
## -(rfe) imp Fold1.Rep2 
## +(rfe) fit Fold1.Rep2 size: 32 
## -(rfe) fit Fold1.Rep2 size: 32 
## +(rfe) fit Fold1.Rep2 size: 16 
## -(rfe) fit Fold1.Rep2 size: 16 
## +(rfe) fit Fold1.Rep2 size:  8 
## -(rfe) fit Fold1.Rep2 size:  8 
## +(rfe) fit Fold1.Rep2 size:  4 
## -(rfe) fit Fold1.Rep2 size:  4 
## +(rfe) fit Fold1.Rep2 size:  2 
## -(rfe) fit Fold1.Rep2 size:  2 
## +(rfe) fit Fold2.Rep2 size: 55
```

```
## Warning in lda.default(x, grouping, ...): variables are collinear
```

```
## -(rfe) fit Fold2.Rep2 size: 55 
## +(rfe) imp Fold2.Rep2 
## -(rfe) imp Fold2.Rep2 
## +(rfe) fit Fold2.Rep2 size: 32 
## -(rfe) fit Fold2.Rep2 size: 32 
## +(rfe) fit Fold2.Rep2 size: 16 
## -(rfe) fit Fold2.Rep2 size: 16 
## +(rfe) fit Fold2.Rep2 size:  8 
## -(rfe) fit Fold2.Rep2 size:  8 
## +(rfe) fit Fold2.Rep2 size:  4 
## -(rfe) fit Fold2.Rep2 size:  4 
## +(rfe) fit Fold2.Rep2 size:  2 
## -(rfe) fit Fold2.Rep2 size:  2 
## +(rfe) fit Fold3.Rep2 size: 55
```

```
## Warning in lda.default(x, grouping, ...): variables are collinear
```

```
## -(rfe) fit Fold3.Rep2 size: 55 
## +(rfe) imp Fold3.Rep2 
## -(rfe) imp Fold3.Rep2 
## +(rfe) fit Fold3.Rep2 size: 32 
## -(rfe) fit Fold3.Rep2 size: 32 
## +(rfe) fit Fold3.Rep2 size: 16 
## -(rfe) fit Fold3.Rep2 size: 16 
## +(rfe) fit Fold3.Rep2 size:  8 
## -(rfe) fit Fold3.Rep2 size:  8 
## +(rfe) fit Fold3.Rep2 size:  4 
## -(rfe) fit Fold3.Rep2 size:  4 
## +(rfe) fit Fold3.Rep2 size:  2 
## -(rfe) fit Fold3.Rep2 size:  2 
## +(rfe) fit Fold1.Rep3 size: 55
```

```
## Warning in lda.default(x, grouping, ...): variables are collinear
```

```
## -(rfe) fit Fold1.Rep3 size: 55 
## +(rfe) imp Fold1.Rep3 
## -(rfe) imp Fold1.Rep3 
## +(rfe) fit Fold1.Rep3 size: 32 
## -(rfe) fit Fold1.Rep3 size: 32 
## +(rfe) fit Fold1.Rep3 size: 16 
## -(rfe) fit Fold1.Rep3 size: 16 
## +(rfe) fit Fold1.Rep3 size:  8 
## -(rfe) fit Fold1.Rep3 size:  8 
## +(rfe) fit Fold1.Rep3 size:  4 
## -(rfe) fit Fold1.Rep3 size:  4 
## +(rfe) fit Fold1.Rep3 size:  2 
## -(rfe) fit Fold1.Rep3 size:  2 
## +(rfe) fit Fold2.Rep3 size: 55
```

```
## Warning in lda.default(x, grouping, ...): variables are collinear
```

```
## -(rfe) fit Fold2.Rep3 size: 55 
## +(rfe) imp Fold2.Rep3 
## -(rfe) imp Fold2.Rep3 
## +(rfe) fit Fold2.Rep3 size: 32 
## -(rfe) fit Fold2.Rep3 size: 32 
## +(rfe) fit Fold2.Rep3 size: 16 
## -(rfe) fit Fold2.Rep3 size: 16 
## +(rfe) fit Fold2.Rep3 size:  8 
## -(rfe) fit Fold2.Rep3 size:  8 
## +(rfe) fit Fold2.Rep3 size:  4 
## -(rfe) fit Fold2.Rep3 size:  4 
## +(rfe) fit Fold2.Rep3 size:  2 
## -(rfe) fit Fold2.Rep3 size:  2 
## +(rfe) fit Fold3.Rep3 size: 55
```

```
## Warning in lda.default(x, grouping, ...): variables are collinear
```

```
## -(rfe) fit Fold3.Rep3 size: 55 
## +(rfe) imp Fold3.Rep3 
## -(rfe) imp Fold3.Rep3 
## +(rfe) fit Fold3.Rep3 size: 32 
## -(rfe) fit Fold3.Rep3 size: 32 
## +(rfe) fit Fold3.Rep3 size: 16 
## -(rfe) fit Fold3.Rep3 size: 16 
## +(rfe) fit Fold3.Rep3 size:  8 
## -(rfe) fit Fold3.Rep3 size:  8 
## +(rfe) fit Fold3.Rep3 size:  4 
## -(rfe) fit Fold3.Rep3 size:  4 
## +(rfe) fit Fold3.Rep3 size:  2 
## -(rfe) fit Fold3.Rep3 size:  2
```

```
## Warning in lda.default(x, grouping, ...): variables are collinear
```

```
## 
## Recursive feature selection
## 
## Outer resampling method: Cross-Validated (3 fold, repeated 3 times) 
## 
## Resampling performance over subset size:
## 
##  Variables Accuracy   Kappa AccuracySD KappaSD Selected
##          2   0.8096 0.03708   0.005112 0.01788         
##          4   0.8864 0.52518   0.005536 0.02625         
##          8   0.8945 0.57850   0.013874 0.06859         
##         16   0.9305 0.75939   0.004671 0.01720         
##         32   0.9303 0.75899   0.004681 0.01717         
##         55   0.9326 0.76882   0.004814 0.01705        *
## 
## The top 5 variables (out of 55):
##    WordCount.log1p, WordCount.root2, WordCount.nexp, NDSS.my.fctrOpEd#Opnn#, PubDate.day.minutes.poly.1
## 
##  [1] "WordCount.log1p"                    
##  [2] "WordCount.root2"                    
##  [3] "WordCount.nexp"                     
##  [4] "NDSS.my.fctrOpEd#Opnn#"             
##  [5] "PubDate.day.minutes.poly.1"         
##  [6] "PubDate.day.minutes.poly.4"         
##  [7] "PubDate.hour.fctr(15.3,23]"         
##  [8] "NDSS.my.fctrScnc#Hlth#"             
##  [9] "PubDate.last4.log1p"                
## [10] "PubDate.last2.log1p"                
## [11] "NDSS.my.fctrBsnss#Crsswrds/Gms#"    
## [12] "NDSS.my.fctrStyls#U.S.#"            
## [13] "PubDate.last8.log1p"                
## [14] "PubDate.day.minutes.poly.5"         
## [15] "PubDate.wkend"                      
## [16] "PubDate.last16.log1p"               
## [17] "PubDate.juliandate"                 
## [18] "PubDate.month.fctr11"               
## [19] "PubDate.day.minutes.poly.3"         
## [20] "PubDate.wkday.fctr6"                
## [21] "PubDate.date.fctr(7,13]"            
## [22] "PubDate.second.fctr(14.8,29.5]"     
## [23] "PubDate.month.fctr10"               
## [24] "PubDate.wkday.fctr1"                
## [25] ".rnorm"                             
## [26] "PubDate.last32.log1p"               
## [27] "PubDate.minute.fctr(44.2,59.1]"     
## [28] "PubDate.day.minutes.poly.2"         
## [29] "PubDate.hour.fctr(7.67,15.3]"       
## [30] "PubDate.minute.fctr(14.8,29.5]"     
## [31] "PubDate.date.fctr(25,31]"           
## [32] "PubDate.second.fctr(44.2,59.1]"     
## [33] "PubDate.wkday.fctr3"                
## [34] "NDSS.my.fctrmyOthr"                 
## [35] "PubDate.date.fctr(19,25]"           
## [36] "NDSS.my.fctr#Opnn#RmFrDbt"          
## [37] "NDSS.my.fctrBsnss#Tchnlgy#"         
## [38] "PubDate.wkday.fctr4"                
## [39] "PubDate.second.fctr(29.5,44.2]"     
## [40] "PubDate.date.fctr(13,19]"           
## [41] "NDSS.my.fctrMtr#N.Y./Rgn#"          
## [42] "NDSS.my.fctrTrvl#Trvl#"             
## [43] "NDSS.my.fctrBsnss#BsnssDy#SmllBsnss"
## [44] "NDSS.my.fctr#Mltmd#"                
## [45] "PubDate.wkday.fctr2"                
## [46] "NDSS.my.fctrStyls##Fshn"            
## [47] "NDSS.my.fctrFrgn#Wrld#"             
## [48] "PubDate.minute.fctr(29.5,44.2]"     
## [49] "NDSS.my.fctrFrgn#Wrld#AsPcfc"       
## [50] "PubDate.wkday.fctr5"                
## [51] "NDSS.my.fctr#U.S.#Edctn"            
## [52] "NDSS.my.fctrCltr#Arts#"             
## [53] "NDSS.my.fctrBsnss#BsnssDy#Dlbk"     
## [54] "NDSS.my.fctr##"                     
## [55] "NDSS.my.fctrTStyl##"
```

![](NYTBlogs3_feat_PubDate_files/figure-html/select.features-3.png) 

```
## [1] "numeric data missing in : "
##   Popular Pplr.fctr 
##      1870      1870 
## [1] "numeric data w/ 0s in : "
##            WordCount              Popular      WordCount.log1p 
##                  109                 5439                  109 
##      WordCount.root2       WordCount.nexp   PubDate.wkday.fctr 
##                  109                 2044                  378 
##        PubDate.wkend        PubDate.hlday  PubDate.day.minutes 
##                 7787                 8160                    5 
##  PubDate.last2.log1p  PubDate.last4.log1p  PubDate.last8.log1p 
##                    2                    4                    8 
## PubDate.last16.log1p PubDate.last32.log1p 
##                   16                   32 
## [1] "numeric data w/ Infs in : "
## named integer(0)
## [1] "numeric data w/ NaNs in : "
## named integer(0)
## [1] "string data missing in : "
##       NewsDesk    SectionName SubsectionName       Headline        Snippet 
##           2408           2899           6176              0             13 
##       Abstract        PubDate        NDSS.my           .lcn 
##             17              0              0           1870
```

```
## [1] "glb_feats_df:"
```

```
## [1] 31 12
```

```
##                  id exclude.as.feat rsp_var
## Pplr.fctr Pplr.fctr            TRUE    TRUE
```

```
##                  id      cor.y exclude.as.feat  cor.y.abs cor.high.X
## Popular     Popular 1.00000000            TRUE 1.00000000       <NA>
## UniqueID   UniqueID 0.01182492            TRUE 0.01182492       <NA>
## Pplr.fctr Pplr.fctr         NA            TRUE         NA       <NA>
##           freqRatio percentUnique zeroVar   nzv is.cor.y.abs.low
## Popular    4.976212    0.03061849   FALSE FALSE            FALSE
## UniqueID   1.000000  100.00000000   FALSE FALSE            FALSE
## Pplr.fctr        NA            NA      NA    NA               NA
##           interaction.feat shapiro.test.p.value rsp_var_raw id_var rsp_var
## Popular                 NA                   NA        TRUE     NA      NA
## UniqueID                NA                   NA       FALSE   TRUE      NA
## Pplr.fctr               NA                   NA          NA     NA    TRUE
```

```
## [1] "glb_feats_df vs. glbObsAll: "
```

```
## character(0)
```

```
## [1] "glbObsAll vs. glb_feats_df: "
```

```
## character(0)
```

```
##              label step_major step_minor label_minor    bgn   end elapsed
## 9  select.features          5          0           0 57.659 80.36  22.701
## 10      fit.models          6          0           0 80.360    NA      NA
```

## Step `6.0: fit models`

```r
fit.models_0_chunk_df <- myadd_chunk(NULL, "fit.models_0_bgn", label.minor = "setup")
```

```
##              label step_major step_minor label_minor    bgn end elapsed
## 1 fit.models_0_bgn          1          0       setup 81.459  NA      NA
```

```r
# load(paste0(glb_out_pfx, "dsk.RData"))

get_model_sel_frmla <- function() {
    model_evl_terms <- c(NULL)
    # min.aic.fit might not be avl
    lclMdlEvlCriteria <- 
        glbMdlMetricsEval[glbMdlMetricsEval %in% names(glb_models_df)]
    for (metric in lclMdlEvlCriteria)
        model_evl_terms <- c(model_evl_terms, 
                             ifelse(length(grep("max", metric)) > 0, "-", "+"), metric)
    if (glb_is_classification && glb_is_binomial)
        model_evl_terms <- c(model_evl_terms, "-", "opt.prob.threshold.OOB")
    model_sel_frmla <- as.formula(paste(c("~ ", model_evl_terms), collapse = " "))
    return(model_sel_frmla)
}

get_dsp_models_df <- function() {
    dsp_models_cols <- c("id", 
                    glbMdlMetricsEval[glbMdlMetricsEval %in% names(glb_models_df)],
                    grep("opt.", names(glb_models_df), fixed = TRUE, value = TRUE)) 
    dsp_models_df <- 
        #orderBy(get_model_sel_frmla(), glb_models_df)[, c("id", glbMdlMetricsEval)]
        orderBy(get_model_sel_frmla(), glb_models_df)[, dsp_models_cols]    
    nCvMdl <- sapply(glb_models_lst, function(mdl) nrow(mdl$results))
    nParams <- sapply(glb_models_lst, function(mdl) ifelse(mdl$method == "custom", 0, 
        nrow(subset(modelLookup(mdl$method), parameter != "parameter"))))
    
#     nCvMdl <- nCvMdl[names(nCvMdl) != "avNNet"]
#     nParams <- nParams[names(nParams) != "avNNet"]    
    
    if (length(cvMdlProblems <- nCvMdl[nCvMdl <= nParams]) > 0) {
        print("Cross Validation issues:")
        warning("Cross Validation issues:")        
        print(cvMdlProblems)
    }
    
    pltMdls <- setdiff(names(nCvMdl), names(cvMdlProblems))
    pltMdls <- setdiff(pltMdls, names(nParams[nParams == 0]))
    
    # length(pltMdls) == 21
    png(paste0(glb_out_pfx, "bestTune.png"), width = 480 * 2, height = 480 * 4)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(ceiling(length(pltMdls) / 2.0), 2)))
    pltIx <- 1
    for (mdlId in pltMdls) {
        print(ggplot(glb_models_lst[[mdlId]], highBestTune = TRUE) + labs(title = mdlId),   
              vp = viewport(layout.pos.row = ceiling(pltIx / 2.0), 
                            layout.pos.col = ((pltIx - 1) %% 2) + 1))  
        pltIx <- pltIx + 1
    }
    dev.off()

    return(dsp_models_df)
}
#get_dsp_models_df()

if (glb_is_classification && glb_is_binomial && 
        (length(unique(glbObsFit[, glb_rsp_var])) < 2))
    stop("glbObsFit$", glb_rsp_var, ": contains less than 2 unique values: ",
         paste0(unique(glbObsFit[, glb_rsp_var]), collapse=", "))

max_cor_y_x_vars <- orderBy(~ -cor.y.abs, 
        subset(glb_feats_df, (exclude.as.feat == 0) & !nzv & !is.cor.y.abs.low & 
                                is.na(cor.high.X)))[1:2, "id"]
max_cor_y_x_vars <- max_cor_y_x_vars[!is.na(max_cor_y_x_vars)]

if (!is.null(glb_Baseline_mdl_var)) {
    if ((max_cor_y_x_vars[1] != glb_Baseline_mdl_var) & 
        (glb_feats_df[glb_feats_df$id == max_cor_y_x_vars[1], "cor.y.abs"] > 
         glb_feats_df[glb_feats_df$id == glb_Baseline_mdl_var, "cor.y.abs"]))
        stop(max_cor_y_x_vars[1], " has a higher correlation with ", glb_rsp_var, 
             " than the Baseline var: ", glb_Baseline_mdl_var)
}

glb_model_type <- ifelse(glb_is_regression, "regression", "classification")
    
# Model specs
c("id.prefix", "method", "type",
  # trainControl params
  "preProc.method", "cv.n.folds", "cv.n.repeats", "summary.fn",
  # train params
  "metric", "metric.maximize", "tune.df")
```

```
##  [1] "id.prefix"       "method"          "type"           
##  [4] "preProc.method"  "cv.n.folds"      "cv.n.repeats"   
##  [7] "summary.fn"      "metric"          "metric.maximize"
## [10] "tune.df"
```

```r
# Baseline
if (!is.null(glb_Baseline_mdl_var)) {
    fit.models_0_chunk_df <- myadd_chunk(fit.models_0_chunk_df, 
                            paste0("fit.models_0_", "Baseline"), major.inc = FALSE,
                                    label.minor = "mybaseln_classfr")
    ret_lst <- myfit_mdl(mdl_id="Baseline", 
                         model_method="mybaseln_classfr",
                        indep_vars_vctr=glb_Baseline_mdl_var,
                        rsp_var=glb_rsp_var,
                        fit_df=glbObsFit, OOB_df=glbObsOOB)
}    

# Most Frequent Outcome "MFO" model: mean(y) for regression
#   Not using caret's nullModel since model stats not avl
#   Cannot use rpart for multinomial classification since it predicts non-MFO
fit.models_0_chunk_df <- myadd_chunk(fit.models_0_chunk_df, 
                            paste0("fit.models_0_", "MFO"), major.inc = FALSE,
                                    label.minor = "myMFO_classfr")
```

```
##              label step_major step_minor   label_minor    bgn    end
## 1 fit.models_0_bgn          1          0         setup 81.459 81.488
## 2 fit.models_0_MFO          1          1 myMFO_classfr 81.489     NA
##   elapsed
## 1    0.03
## 2      NA
```

```r
ret_lst <- myfit_mdl(mdl_specs_lst = myinit_mdl_specs_lst(mdl_specs_lst = list(
    id.prefix = "MFO", type = glb_model_type, trainControl.method = "none",
    train.method = ifelse(glb_is_regression, "lm", "myMFO_classfr"))),
                        indep_vars = ".rnorm", rsp_var = glb_rsp_var,
                        fit_df = glbObsFit, OOB_df = glbObsOOB)
```

```
## [1] "fitting model: MFO###myMFO_classfr"
## [1] "    indep_vars: .rnorm"
## Fitting parameter = none on full training set
## [1] "in MFO.Classifier$fit"
## [1] "unique.vals:"
## [1] N Y
## Levels: N Y
## [1] "unique.prob:"
## y
##        N        Y 
## 0.820358 0.179642 
## [1] "MFO.val:"
## [1] "N"
##             Length Class      Mode     
## unique.vals 2      factor     numeric  
## unique.prob 2      -none-     numeric  
## MFO.val     1      -none-     character
## x.names     1      -none-     character
## xNames      1      -none-     character
## problemType 1      -none-     character
## tuneValue   1      data.frame list     
## obsLevels   2      -none-     character
## [1] "entr MFO.Classifier$predict"
## [1] "exit MFO.Classifier$predict"
```

```
## Loading required package: ROCR
## Loading required package: gplots
## 
## Attaching package: 'gplots'
## 
## The following object is masked from 'package:stats':
## 
##     lowess
```

```
## [1] "in MFO.Classifier$prob"
##          N        Y
## 1 0.820358 0.179642
## 2 0.820358 0.179642
## 3 0.820358 0.179642
## 4 0.820358 0.179642
## 5 0.820358 0.179642
## 6 0.820358 0.179642
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-1.png) 

```
##          Prediction
## Reference    N    Y
##         N    0 3941
##         Y    0  863
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##      0.1796420      0.0000000      0.1688795      0.1907952      0.8203580 
## AccuracyPValue  McnemarPValue 
##      1.0000000      0.0000000 
## [1] "entr MFO.Classifier$predict"
## [1] "exit MFO.Classifier$predict"
## [1] "in MFO.Classifier$prob"
##          N        Y
## 1 0.820358 0.179642
## 2 0.820358 0.179642
## 3 0.820358 0.179642
## 4 0.820358 0.179642
## 5 0.820358 0.179642
## 6 0.820358 0.179642
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-2.png) 

```
##          Prediction
## Reference    N    Y
##         N    0 1498
##         Y    0  230
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##      0.1331019      0.0000000      0.1174298      0.1500310      0.8668981 
## AccuracyPValue  McnemarPValue 
##      1.0000000      0.0000000 
##                    id  feats max.nTuningRuns min.elapsedtime.everything
## 1 MFO###myMFO_classfr .rnorm               0                      0.302
##   min.elapsedtime.final max.AUCpROC.fit max.Sens.fit max.Spec.fit
## 1                 0.004             0.5            1            0
##   max.AUCROCR.fit opt.prob.threshold.fit max.f.score.fit max.Accuracy.fit
## 1             0.5                    0.1       0.3045703         0.179642
##   max.AccuracyLower.fit max.AccuracyUpper.fit max.Kappa.fit
## 1             0.1688795             0.1907952             0
##   max.AUCpROC.OOB max.Sens.OOB max.Spec.OOB max.AUCROCR.OOB
## 1             0.5            1            0             0.5
##   opt.prob.threshold.OOB max.f.score.OOB max.Accuracy.OOB
## 1                    0.1       0.2349336        0.1331019
##   max.AccuracyLower.OOB max.AccuracyUpper.OOB max.Kappa.OOB
## 1             0.1174298              0.150031             0
```

```r
if (glb_is_classification) {
    # "random" model - only for classification; 
    #   none needed for regression since it is same as MFO
fit.models_0_chunk_df <- myadd_chunk(fit.models_0_chunk_df, 
                            paste0("fit.models_0_", "Random"), major.inc = FALSE,
                                    label.minor = "myrandom_classfr")

#stop(here"); glb2Sav(); all.equal(glb_models_df, sav_models_df)    
    ret_lst <- myfit_mdl(mdl_specs_lst = myinit_mdl_specs_lst(mdl_specs_lst = list(
        id.prefix = "Random", type = glb_model_type, trainControl.method = "none",
        train.method = "myrandom_classfr")),
                        indep_vars = ".rnorm", rsp_var = glb_rsp_var,
                        fit_df = glbObsFit, OOB_df = glbObsOOB)
}
```

```
##                 label step_major step_minor      label_minor    bgn    end
## 2    fit.models_0_MFO          1          1    myMFO_classfr 81.489 84.531
## 3 fit.models_0_Random          1          2 myrandom_classfr 84.532     NA
##   elapsed
## 2   3.042
## 3      NA
## [1] "fitting model: Random###myrandom_classfr"
## [1] "    indep_vars: .rnorm"
## Fitting parameter = none on full training set
##             Length Class      Mode     
## unique.vals 2      factor     numeric  
## unique.prob 2      table      numeric  
## xNames      1      -none-     character
## problemType 1      -none-     character
## tuneValue   1      data.frame list     
## obsLevels   2      -none-     character
## [1] "in Random.Classifier$prob"
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-3.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-4.png) 

```
##          Prediction
## Reference    N    Y
##         N    0 3941
##         Y    0  863
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##      0.1796420      0.0000000      0.1688795      0.1907952      0.8203580 
## AccuracyPValue  McnemarPValue 
##      1.0000000      0.0000000 
## [1] "in Random.Classifier$prob"
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-5.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-6.png) 

```
##          Prediction
## Reference    N    Y
##         N    0 1498
##         Y    0  230
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##      0.1331019      0.0000000      0.1174298      0.1500310      0.8668981 
## AccuracyPValue  McnemarPValue 
##      1.0000000      0.0000000 
##                          id  feats max.nTuningRuns
## 1 Random###myrandom_classfr .rnorm               0
##   min.elapsedtime.everything min.elapsedtime.final max.AUCpROC.fit
## 1                      0.309                 0.002       0.4990604
##   max.Sens.fit max.Spec.fit max.AUCROCR.fit opt.prob.threshold.fit
## 1    0.8312611    0.1668598       0.4972757                    0.1
##   max.f.score.fit max.Accuracy.fit max.AccuracyLower.fit
## 1       0.3045703         0.179642             0.1688795
##   max.AccuracyUpper.fit max.Kappa.fit max.AUCpROC.OOB max.Sens.OOB
## 1             0.1907952             0       0.5125675    0.8077437
##   max.Spec.OOB max.AUCROCR.OOB opt.prob.threshold.OOB max.f.score.OOB
## 1    0.2173913       0.4857956                    0.1       0.2349336
##   max.Accuracy.OOB max.AccuracyLower.OOB max.AccuracyUpper.OOB
## 1        0.1331019             0.1174298              0.150031
##   max.Kappa.OOB
## 1             0
```

```r
# Max.cor.Y
#   Check impact of cv
#       rpart is not a good candidate since caret does not optimize cp (only tuning parameter of rpart) well
fit.models_0_chunk_df <- myadd_chunk(fit.models_0_chunk_df, 
                        paste0("fit.models_0_", "Max.cor.Y.rcv.*X*"), major.inc = FALSE,
                                    label.minor = "glmnet")
```

```
##                            label step_major step_minor      label_minor
## 3            fit.models_0_Random          1          2 myrandom_classfr
## 4 fit.models_0_Max.cor.Y.rcv.*X*          1          3           glmnet
##      bgn    end elapsed
## 3 84.532 88.901    4.37
## 4 88.902     NA      NA
```

```r
ret_lst <- myfit_mdl(mdl_specs_lst=myinit_mdl_specs_lst(mdl_specs_lst=list(
    id.prefix="Max.cor.Y.rcv.1X1", type=glb_model_type, trainControl.method="none",
    train.method="glmnet")),
                    indep_vars=max_cor_y_x_vars, rsp_var=glb_rsp_var, 
                    fit_df=glbObsFit, OOB_df=glbObsOOB)
```

```
## [1] "fitting model: Max.cor.Y.rcv.1X1###glmnet"
## [1] "    indep_vars: WordCount.root2,NDSS.my.fctr"
```

```
## Loading required package: glmnet
## Loading required package: Matrix
## Loaded glmnet 2.0-2
```

```
## Fitting alpha = 0.1, lambda = 0.00434 on full training set
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-7.png) 

```
##             Length Class      Mode     
## a0           100   -none-     numeric  
## beta        2100   dgCMatrix  S4       
## df           100   -none-     numeric  
## dim            2   -none-     numeric  
## lambda       100   -none-     numeric  
## dev.ratio    100   -none-     numeric  
## nulldev        1   -none-     numeric  
## npasses        1   -none-     numeric  
## jerr           1   -none-     numeric  
## offset         1   -none-     logical  
## classnames     2   -none-     character
## call           5   -none-     call     
## nobs           1   -none-     numeric  
## lambdaOpt      1   -none-     numeric  
## xNames        21   -none-     character
## problemType    1   -none-     character
## tuneValue      2   data.frame list     
## obsLevels      2   -none-     character
## [1] "min lambda > lambdaOpt:"
##                         (Intercept)                 NDSS.my.fctr#Mltmd# 
##                         -4.57159198                         -1.22219085 
##           NDSS.my.fctr#Opnn#RmFrDbt        NDSS.my.fctr#Opnn#ThPblcEdtr 
##                         -3.46072453                          4.06871185 
##             NDSS.my.fctr#U.S.#Edctn      NDSS.my.fctrBsnss#BsnssDy#Dlbk 
##                         -1.89443632                         -0.22472818 
## NDSS.my.fctrBsnss#BsnssDy#SmllBsnss     NDSS.my.fctrBsnss#Crsswrds/Gms# 
##                         -0.95537118                          4.55408513 
##          NDSS.my.fctrBsnss#Tchnlgy#              NDSS.my.fctrCltr#Arts# 
##                          0.77368538                         -0.09465691 
##              NDSS.my.fctrFrgn#Wrld#        NDSS.my.fctrFrgn#Wrld#AsPcfc 
##                         -1.45528874                         -1.60117505 
##           NDSS.my.fctrMtr#N.Y./Rgn#              NDSS.my.fctrOpEd#Opnn# 
##                          0.01563989                          4.51696382 
##              NDSS.my.fctrScnc#Hlth#             NDSS.my.fctrStyls##Fshn 
##                          3.51595317                         -1.85948925 
##             NDSS.my.fctrStyls#U.S.#                 NDSS.my.fctrTStyl## 
##                          3.27995325                         -1.54110404 
##              NDSS.my.fctrTrvl#Trvl#                  NDSS.my.fctrmyOthr 
##                         -1.41940605                         -1.90156922 
##                     WordCount.root2 
##                          0.08434378 
## [1] "max lambda < lambdaOpt:"
##                         (Intercept)                 NDSS.my.fctr#Mltmd# 
##                         -4.60394059                         -1.25163328 
##           NDSS.my.fctr#Opnn#RmFrDbt        NDSS.my.fctr#Opnn#ThPblcEdtr 
##                         -3.55521332                          4.09217313 
##             NDSS.my.fctr#U.S.#Edctn      NDSS.my.fctrBsnss#BsnssDy#Dlbk 
##                         -1.96172971                         -0.22495986 
## NDSS.my.fctrBsnss#BsnssDy#SmllBsnss     NDSS.my.fctrBsnss#Crsswrds/Gms# 
##                         -0.96836050                          4.58120497 
##          NDSS.my.fctrBsnss#Tchnlgy#              NDSS.my.fctrCltr#Arts# 
##                          0.78504703                         -0.09069661 
##              NDSS.my.fctrFrgn#Wrld#        NDSS.my.fctrFrgn#Wrld#AsPcfc 
##                         -1.51061232                         -1.63313235 
##           NDSS.my.fctrMtr#N.Y./Rgn#              NDSS.my.fctrOpEd#Opnn# 
##                          0.02466697                          4.54361134 
##              NDSS.my.fctrScnc#Hlth#             NDSS.my.fctrStyls##Fshn 
##                          3.53210055                         -1.92188290 
##             NDSS.my.fctrStyls#U.S.#                 NDSS.my.fctrTStyl## 
##                          3.29488750                         -1.57788931 
##              NDSS.my.fctrTrvl#Trvl#                  NDSS.my.fctrmyOthr 
##                         -1.47368131                         -1.97357582 
##                     WordCount.root2 
##                          0.08537319
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-8.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-9.png) 

```
##          Prediction
## Reference    N    Y
##         N 3796  145
##         Y  177  686
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   9.329725e-01   7.692476e-01   9.255302e-01   9.398832e-01   8.203580e-01 
## AccuracyPValue  McnemarPValue 
##  1.026390e-114   8.406670e-02
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-10.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-11.png) 

```
##          Prediction
## Reference    N    Y
##         N 1151  347
##         Y   67  163
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   7.604167e-01   3.148374e-01   7.395703e-01   7.803749e-01   8.668981e-01 
## AccuracyPValue  McnemarPValue 
##   1.000000e+00   8.593187e-43 
##                           id                        feats max.nTuningRuns
## 1 Max.cor.Y.rcv.1X1###glmnet WordCount.root2,NDSS.my.fctr               0
##   min.elapsedtime.everything min.elapsedtime.final max.AUCpROC.fit
## 1                      1.034                 0.282       0.8790544
##   max.Sens.fit max.Spec.fit max.AUCROCR.fit opt.prob.threshold.fit
## 1    0.9632073    0.7949015       0.9608594                    0.5
##   max.f.score.fit max.Accuracy.fit max.AccuracyLower.fit
## 1       0.8099174        0.9329725             0.9255302
##   max.AccuracyUpper.fit max.Kappa.fit max.AUCpROC.OOB max.Sens.OOB
## 1             0.9398832     0.7692476       0.5962443    0.9098798
##   max.Spec.OOB max.AUCROCR.OOB opt.prob.threshold.OOB max.f.score.OOB
## 1    0.2826087       0.8116126                    0.1       0.4405405
##   max.Accuracy.OOB max.AccuracyLower.OOB max.AccuracyUpper.OOB
## 1        0.7604167             0.7395703             0.7803749
##   max.Kappa.OOB
## 1     0.3148374
```

```r
if (glbMdlCheckRcv) {
    # rcv_n_folds == 1 & rcv_n_repeats > 1 crashes
    for (rcv_n_folds in seq(3, glb_rcv_n_folds + 2, 2))
        for (rcv_n_repeats in seq(1, glb_rcv_n_repeats + 2, 2)) {
            
            # Experiment specific code to avoid caret crash
    #         lcl_tune_models_df <- rbind(data.frame()
    #                             ,data.frame(method = "glmnet", parameter = "alpha", 
    #                                         vals = "0.100 0.325 0.550 0.775 1.000")
    #                             ,data.frame(method = "glmnet", parameter = "lambda",
    #                                         vals = "9.342e-02")    
    #                                     )
            
            ret_lst <- myfit_mdl(mdl_specs_lst = myinit_mdl_specs_lst(mdl_specs_lst =
                list(
                id.prefix = paste0("Max.cor.Y.rcv.", rcv_n_folds, "X", rcv_n_repeats), 
                type = glb_model_type, 
    # tune.df = lcl_tune_models_df,            
                trainControl.method = "repeatedcv",
                trainControl.number = rcv_n_folds, 
                trainControl.repeats = rcv_n_repeats,
                trainControl.classProbs = glb_is_classification,
                trainControl.summaryFunction = glbMdlMetricSummaryFn,
                train.method = "glmnet", train.metric = glbMdlMetricSummary, 
                train.maximize = glbMdlMetricMaximize)),
                                indep_vars = max_cor_y_x_vars, rsp_var = glb_rsp_var, 
                                fit_df = glbObsFit, OOB_df = glbObsOOB)
        }
    # Add parallel coordinates graph of glb_models_df[, glbMdlMetricsEval] to evaluate cv parameters
    tmp_models_cols <- c("id", "max.nTuningRuns",
                        glbMdlMetricsEval[glbMdlMetricsEval %in% names(glb_models_df)],
                        grep("opt.", names(glb_models_df), fixed = TRUE, value = TRUE)) 
    print(myplot_parcoord(obs_df = subset(glb_models_df, 
                                          grepl("Max.cor.Y.rcv.", id, fixed = TRUE), 
                                            select = -feats)[, tmp_models_cols],
                          id_var = "id"))
}

# Useful for stacking decisions
# fit.models_0_chunk_df <- myadd_chunk(fit.models_0_chunk_df, 
#                     paste0("fit.models_0_", "Max.cor.Y[rcv.1X1.cp.0|]"), major.inc = FALSE,
#                                     label.minor = "rpart")
# 
# ret_lst <- myfit_mdl(mdl_specs_lst = myinit_mdl_specs_lst(mdl_specs_lst = list(
#     id.prefix = "Max.cor.Y.rcv.1X1.cp.0", type = glb_model_type, trainControl.method = "none",
#     train.method = "rpart",
#     tune.df=data.frame(method="rpart", parameter="cp", min=0.0, max=0.0, by=0.1))),
#                     indep_vars=max_cor_y_x_vars, rsp_var=glb_rsp_var, 
#                     fit_df=glbObsFit, OOB_df=glbObsOOB)

#stop(here"); glb2Sav(); all.equal(glb_models_df, sav_models_df)
# if (glb_is_regression || glb_is_binomial) # For multinomials this model will be run next by default
ret_lst <- myfit_mdl(mdl_specs_lst=myinit_mdl_specs_lst(mdl_specs_lst=list(
    id.prefix="Max.cor.Y", 
    type=glb_model_type, trainControl.method="repeatedcv",
    trainControl.number=glb_rcv_n_folds, trainControl.repeats=glb_rcv_n_repeats,
            trainControl.classProbs = glb_is_classification,
            trainControl.summaryFunction = glbMdlMetricSummaryFn,
trainControl.allowParallel = FALSE,
            train.metric = glbMdlMetricSummary, 
            train.maximize = glbMdlMetricMaximize,    
    train.method="rpart")),
    indep_vars=max_cor_y_x_vars, rsp_var=glb_rsp_var, 
    fit_df=glbObsFit, OOB_df=glbObsOOB)
```

```
## [1] "fitting model: Max.cor.Y##rcv#rpart"
## [1] "    indep_vars: WordCount.root2,NDSS.my.fctr"
```

```
## Loading required package: rpart
```

```
## + Fold1.Rep1: cp=0.01043 
## - Fold1.Rep1: cp=0.01043 
## + Fold2.Rep1: cp=0.01043 
## - Fold2.Rep1: cp=0.01043 
## + Fold3.Rep1: cp=0.01043 
## - Fold3.Rep1: cp=0.01043 
## + Fold1.Rep2: cp=0.01043 
## - Fold1.Rep2: cp=0.01043 
## + Fold2.Rep2: cp=0.01043 
## - Fold2.Rep2: cp=0.01043 
## + Fold3.Rep2: cp=0.01043 
## - Fold3.Rep2: cp=0.01043 
## + Fold1.Rep3: cp=0.01043 
## - Fold1.Rep3: cp=0.01043 
## + Fold2.Rep3: cp=0.01043 
## - Fold2.Rep3: cp=0.01043 
## + Fold3.Rep3: cp=0.01043 
## - Fold3.Rep3: cp=0.01043 
## Aggregating results
## Selecting tuning parameters
## Fitting cp = 0.0104 on full training set
```

```
## Warning in myfit_mdl(mdl_specs_lst = myinit_mdl_specs_lst(mdl_specs_lst
## = list(id.prefix = "Max.cor.Y", : model's bestTune found at an extreme of
## tuneGrid for parameter: cp
```

```
## Loading required package: rpart.plot
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-12.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-13.png) 

```
## Call:
## rpart(formula = .outcome ~ ., control = list(minsplit = 20, minbucket = 7, 
##     cp = 0, maxcompete = 4, maxsurrogate = 5, usesurrogate = 2, 
##     surrogatestyle = 0, maxdepth = 30, xval = 0))
##   n= 4804 
## 
##           CP nsplit rel error
## 1 0.36964079      0 1.0000000
## 2 0.09849363      1 0.6303592
## 3 0.08574739      2 0.5318656
## 4 0.05677868      3 0.4461182
## 5 0.01042874      4 0.3893395
## 
## Variable importance
##          NDSS.my.fctrOpEd#Opnn# NDSS.my.fctrBsnss#Crsswrds/Gms# 
##                              55                              16 
##          NDSS.my.fctrScnc#Hlth#         NDSS.my.fctrStyls#U.S.# 
##                              16                              12 
## 
## Node number 1: 4804 observations,    complexity param=0.3696408
##   predicted class=N  expected loss=0.179642  P(node) =1
##     class counts:  3941   863
##    probabilities: 0.820 0.180 
##   left son=2 (4367 obs) right son=3 (437 obs)
##   Primary splits:
##       NDSS.my.fctrOpEd#Opnn#          < 0.5      to the left,  improve=451.59770, (0 missing)
##       NDSS.my.fctrBsnss#Crsswrds/Gms# < 0.5      to the left,  improve=112.88510, (0 missing)
##       WordCount.root2                 < 25.75849 to the left,  improve=111.17610, (0 missing)
##       NDSS.my.fctrScnc#Hlth#          < 0.5      to the left,  improve= 99.35206, (0 missing)
##       NDSS.my.fctrStyls#U.S.#         < 0.5      to the left,  improve= 68.73272, (0 missing)
## 
## Node number 2: 4367 observations,    complexity param=0.09849363
##   predicted class=N  expected loss=0.1110602  P(node) =0.9090341
##     class counts:  3882   485
##    probabilities: 0.889 0.111 
##   left son=4 (4262 obs) right son=5 (105 obs)
##   Primary splits:
##       NDSS.my.fctrBsnss#Crsswrds/Gms# < 0.5      to the left,  improve=135.55130, (0 missing)
##       NDSS.my.fctrScnc#Hlth#          < 0.5      to the left,  improve=125.07920, (0 missing)
##       WordCount.root2                 < 25.75849 to the left,  improve= 94.70710, (0 missing)
##       NDSS.my.fctrStyls#U.S.#         < 0.5      to the left,  improve= 88.56821, (0 missing)
##       NDSS.my.fctr#Opnn#ThPblcEdtr    < 0.5      to the left,  improve= 18.74400, (0 missing)
## 
## Node number 3: 437 observations
##   predicted class=Y  expected loss=0.1350114  P(node) =0.09096586
##     class counts:    59   378
##    probabilities: 0.135 0.865 
## 
## Node number 4: 4262 observations,    complexity param=0.08574739
##   predicted class=N  expected loss=0.09150634  P(node) =0.8871774
##     class counts:  3872   390
##    probabilities: 0.908 0.092 
##   left son=8 (4114 obs) right son=9 (148 obs)
##   Primary splits:
##       NDSS.my.fctrScnc#Hlth#       < 0.5      to the left,  improve=132.96710, (0 missing)
##       NDSS.my.fctrStyls#U.S.#      < 0.5      to the left,  improve= 94.69099, (0 missing)
##       WordCount.root2              < 26.49528 to the left,  improve= 84.07487, (0 missing)
##       NDSS.my.fctr#Opnn#ThPblcEdtr < 0.5      to the left,  improve= 19.71762, (0 missing)
##       NDSS.my.fctrTStyl##          < 0.5      to the right, improve= 10.17000, (0 missing)
## 
## Node number 5: 105 observations
##   predicted class=Y  expected loss=0.0952381  P(node) =0.02185679
##     class counts:    10    95
##    probabilities: 0.095 0.905 
## 
## Node number 8: 4114 observations,    complexity param=0.05677868
##   predicted class=N  expected loss=0.06781721  P(node) =0.8563697
##     class counts:  3835   279
##    probabilities: 0.932 0.068 
##   left son=16 (3987 obs) right son=17 (127 obs)
##   Primary splits:
##       NDSS.my.fctrStyls#U.S.#      < 0.5      to the left,  improve=102.410700, (0 missing)
##       WordCount.root2              < 25.01    to the left,  improve= 47.352210, (0 missing)
##       NDSS.my.fctr#Opnn#ThPblcEdtr < 0.5      to the left,  improve= 20.930810, (0 missing)
##       NDSS.my.fctrTStyl##          < 0.5      to the right, improve=  5.249425, (0 missing)
##       NDSS.my.fctrBsnss#Tchnlgy#   < 0.5      to the left,  improve=  2.395935, (0 missing)
## 
## Node number 9: 148 observations
##   predicted class=Y  expected loss=0.25  P(node) =0.03080766
##     class counts:    37   111
##    probabilities: 0.250 0.750 
## 
## Node number 16: 3987 observations
##   predicted class=N  expected loss=0.04790569  P(node) =0.8299334
##     class counts:  3796   191
##    probabilities: 0.952 0.048 
## 
## Node number 17: 127 observations
##   predicted class=Y  expected loss=0.3070866  P(node) =0.0264363
##     class counts:    39    88
##    probabilities: 0.307 0.693 
## 
## n= 4804 
## 
## node), split, n, loss, yval, (yprob)
##       * denotes terminal node
## 
##  1) root 4804 863 N (0.82035803 0.17964197)  
##    2) NDSS.my.fctrOpEd#Opnn#< 0.5 4367 485 N (0.88893978 0.11106022)  
##      4) NDSS.my.fctrBsnss#Crsswrds/Gms#< 0.5 4262 390 N (0.90849366 0.09150634)  
##        8) NDSS.my.fctrScnc#Hlth#< 0.5 4114 279 N (0.93218279 0.06781721)  
##         16) NDSS.my.fctrStyls#U.S.#< 0.5 3987 191 N (0.95209431 0.04790569) *
##         17) NDSS.my.fctrStyls#U.S.#>=0.5 127  39 Y (0.30708661 0.69291339) *
##        9) NDSS.my.fctrScnc#Hlth#>=0.5 148  37 Y (0.25000000 0.75000000) *
##      5) NDSS.my.fctrBsnss#Crsswrds/Gms#>=0.5 105  10 Y (0.09523810 0.90476190) *
##    3) NDSS.my.fctrOpEd#Opnn#>=0.5 437  59 Y (0.13501144 0.86498856) *
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-14.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-15.png) 

```
##          Prediction
## Reference    N    Y
##         N 3796  145
##         Y  191  672
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   9.300583e-01   7.576571e-01   9.224771e-01   9.371115e-01   8.203580e-01 
## AccuracyPValue  McnemarPValue 
##  4.458834e-108   1.409037e-02
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-16.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-17.png) 

```
##          Prediction
## Reference    N    Y
##         N 1355  143
##         Y  168   62
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##      0.8200231      0.1825002      0.8010821      0.8378705      0.8668981 
## AccuracyPValue  McnemarPValue 
##      1.0000000      0.1735405 
##                     id                        feats max.nTuningRuns
## 1 Max.cor.Y##rcv#rpart WordCount.root2,NDSS.my.fctr               5
##   min.elapsedtime.everything min.elapsedtime.final max.AUCpROC.fit
## 1                      3.029                 0.077       0.8709432
##   max.Sens.fit max.Spec.fit max.AUCROCR.fit opt.prob.threshold.fit
## 1    0.9632073     0.778679       0.8746354                    0.6
##   max.f.score.fit max.Accuracy.fit max.AccuracyLower.fit
## 1             0.8        0.9296422             0.9224771
##   max.AccuracyUpper.fit max.Kappa.fit max.AUCpROC.OOB max.Sens.OOB
## 1             0.9371115     0.7515134       0.5870523    0.9045394
##   max.Spec.OOB max.AUCROCR.OOB opt.prob.threshold.OOB max.f.score.OOB
## 1    0.2695652       0.5892132                    0.6       0.2850575
##   max.Accuracy.OOB max.AccuracyLower.OOB max.AccuracyUpper.OOB
## 1        0.8200231             0.8010821             0.8378705
##   max.Kappa.OOB max.AccuracySD.fit max.KappaSD.fit
## 1     0.1825002         0.00506952       0.0191091
```

```r
if ((length(glbFeatsDateTime) > 0) && 
    (sum(grepl(paste(names(glbFeatsDateTime), "\\.day\\.minutes\\.poly\\.", sep = ""),
               names(glbObsAll))) > 0)) {
    fit.models_0_chunk_df <- myadd_chunk(fit.models_0_chunk_df, 
                    paste0("fit.models_0_", "Max.cor.Y.Time.Poly"), major.inc = FALSE,
                                    label.minor = "glmnet")

    indepVars <- c(max_cor_y_x_vars, 
            grep(paste(names(glbFeatsDateTime), "\\.day\\.minutes\\.poly\\.", sep = ""),
                        names(glbObsAll), value = TRUE))
    indepVars <- myadjust_interaction_feats(indepVars)
    ret_lst <- myfit_mdl(mdl_specs_lst = myinit_mdl_specs_lst(mdl_specs_lst = list(
        id.prefix = "Max.cor.Y.Time.Poly", 
        type = glb_model_type, trainControl.method = "repeatedcv",
        trainControl.number = glb_rcv_n_folds, trainControl.repeats = glb_rcv_n_repeats,
        trainControl.classProbs = glb_is_classification,
        trainControl.summaryFunction = glbMdlMetricSummaryFn,
        train.metric = glbMdlMetricSummary, 
        train.maximize = glbMdlMetricMaximize,    
        train.method = "glmnet")),
        indep_vars = indepVars,
        rsp_var = glb_rsp_var, 
        fit_df = glbObsFit, OOB_df = glbObsOOB)
}
```

```
##                              label step_major step_minor label_minor
## 4   fit.models_0_Max.cor.Y.rcv.*X*          1          3      glmnet
## 5 fit.models_0_Max.cor.Y.Time.Poly          1          4      glmnet
##       bgn     end elapsed
## 4  88.902 102.249  13.347
## 5 102.250      NA      NA
## [1] "fitting model: Max.cor.Y.Time.Poly##rcv#glmnet"
## [1] "    indep_vars: WordCount.root2,NDSS.my.fctr,PubDate.day.minutes.poly.1,PubDate.day.minutes.poly.2,PubDate.day.minutes.poly.3,PubDate.day.minutes.poly.4,PubDate.day.minutes.poly.5"
## Aggregating results
## Selecting tuning parameters
## Fitting alpha = 0.55, lambda = 0.0201 on full training set
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-18.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-19.png) 

```
##             Length Class      Mode     
## a0            95   -none-     numeric  
## beta        2470   dgCMatrix  S4       
## df            95   -none-     numeric  
## dim            2   -none-     numeric  
## lambda        95   -none-     numeric  
## dev.ratio     95   -none-     numeric  
## nulldev        1   -none-     numeric  
## npasses        1   -none-     numeric  
## jerr           1   -none-     numeric  
## offset         1   -none-     logical  
## classnames     2   -none-     character
## call           5   -none-     call     
## nobs           1   -none-     numeric  
## lambdaOpt      1   -none-     numeric  
## xNames        26   -none-     character
## problemType    1   -none-     character
## tuneValue      2   data.frame list     
## obsLevels      2   -none-     character
## [1] "min lambda > lambdaOpt:"
##                     (Intercept)       NDSS.my.fctr#Opnn#RmFrDbt 
##                    -3.824918286                    -0.730902967 
##    NDSS.my.fctr#Opnn#ThPblcEdtr NDSS.my.fctrBsnss#Crsswrds/Gms# 
##                     2.927432404                     3.640472806 
##      NDSS.my.fctrBsnss#Tchnlgy#    NDSS.my.fctrFrgn#Wrld#AsPcfc 
##                     0.192915649                    -0.002462006 
##          NDSS.my.fctrOpEd#Opnn#          NDSS.my.fctrScnc#Hlth# 
##                     3.940210915                     3.120520277 
##         NDSS.my.fctrStyls#U.S.#             NDSS.my.fctrTStyl## 
##                     2.887354523                    -0.349210176 
##      PubDate.day.minutes.poly.1      PubDate.day.minutes.poly.2 
##                    10.171710001                     1.938052563 
##      PubDate.day.minutes.poly.4                 WordCount.root2 
##                     0.422263612                     0.053515500 
## [1] "max lambda < lambdaOpt:"
##                     (Intercept)       NDSS.my.fctr#Opnn#RmFrDbt 
##                    -3.902999004                    -0.872538727 
##    NDSS.my.fctr#Opnn#ThPblcEdtr         NDSS.my.fctr#U.S.#Edctn 
##                     3.038605329                    -0.004334849 
## NDSS.my.fctrBsnss#Crsswrds/Gms#      NDSS.my.fctrBsnss#Tchnlgy# 
##                     3.700053147                     0.274702193 
##    NDSS.my.fctrFrgn#Wrld#AsPcfc          NDSS.my.fctrOpEd#Opnn# 
##                    -0.061833760                     4.010945376 
##          NDSS.my.fctrScnc#Hlth#         NDSS.my.fctrStyls#U.S.# 
##                     3.177975989                     2.950583177 
##             NDSS.my.fctrTStyl##      PubDate.day.minutes.poly.1 
##                    -0.393389148                    11.082450379 
##      PubDate.day.minutes.poly.2      PubDate.day.minutes.poly.4 
##                     2.849459283                     1.090937007 
##                 WordCount.root2 
##                     0.055710821
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-20.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-21.png) 

```
##          Prediction
## Reference    N    Y
##         N 3796  145
##         Y  177  686
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   9.329725e-01   7.692476e-01   9.255302e-01   9.398832e-01   8.203580e-01 
## AccuracyPValue  McnemarPValue 
##  1.026390e-114   8.406670e-02
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-22.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-23.png) 

```
##          Prediction
## Reference    N    Y
##         N 1185  313
##         Y   75  155
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   7.754630e-01   3.233542e-01   7.550404e-01   7.949457e-01   8.668981e-01 
## AccuracyPValue  McnemarPValue 
##   1.000000e+00   2.416777e-33 
##                                id
## 1 Max.cor.Y.Time.Poly##rcv#glmnet
##                                                                                                                                                                 feats
## 1 WordCount.root2,NDSS.my.fctr,PubDate.day.minutes.poly.1,PubDate.day.minutes.poly.2,PubDate.day.minutes.poly.3,PubDate.day.minutes.poly.4,PubDate.day.minutes.poly.5
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1              25                      4.889                 0.316
##   max.AUCpROC.fit max.Sens.fit max.Spec.fit max.AUCROCR.fit
## 1        0.874855    0.9652372    0.7844728       0.9565422
##   opt.prob.threshold.fit max.f.score.fit max.Accuracy.fit
## 1                    0.4       0.8099174         0.932279
##   max.AccuracyLower.fit max.AccuracyUpper.fit max.Kappa.fit
## 1             0.9255302             0.9398832     0.7643948
##   max.AUCpROC.OOB max.Sens.OOB max.Spec.OOB max.AUCROCR.OOB
## 1        0.596578    0.9105474    0.2826087       0.8049472
##   opt.prob.threshold.OOB max.f.score.OOB max.Accuracy.OOB
## 1                    0.1       0.4441261         0.775463
##   max.AccuracyLower.OOB max.AccuracyUpper.OOB max.Kappa.OOB
## 1             0.7550404             0.7949457     0.3233542
##   max.AccuracySD.fit max.KappaSD.fit
## 1         0.00531278      0.01857466
```

```r
if ((length(glbFeatsDateTime) > 0) && 
    (sum(grepl(paste(names(glbFeatsDateTime), "\\.last[[:digit:]]", sep = ""),
               names(glbObsAll))) > 0)) {
    fit.models_0_chunk_df <- myadd_chunk(fit.models_0_chunk_df, 
                    paste0("fit.models_0_", "Max.cor.Y.Time.Lag"), major.inc = FALSE,
                                    label.minor = "glmnet")

    indepVars <- c(max_cor_y_x_vars, 
            grep(paste(names(glbFeatsDateTime), "\\.last[[:digit:]]", sep = ""),
                        names(glbObsAll), value = TRUE))
    indepVars <- myadjust_interaction_feats(indepVars)
    ret_lst <- myfit_mdl(mdl_specs_lst = myinit_mdl_specs_lst(mdl_specs_lst = list(
        id.prefix = "Max.cor.Y.Time.Lag", 
        type = glb_model_type, 
        tune.df = glbMdlTuneParams,        
        trainControl.method = "repeatedcv",
        trainControl.number = glb_rcv_n_folds, trainControl.repeats = glb_rcv_n_repeats,
        trainControl.classProbs = glb_is_classification,
        trainControl.summaryFunction = glbMdlMetricSummaryFn,
        train.metric = glbMdlMetricSummary, 
        train.maximize = glbMdlMetricMaximize,    
        train.method = "glmnet")),
        indep_vars = indepVars,
        rsp_var = glb_rsp_var, 
        fit_df = glbObsFit, OOB_df = glbObsOOB)
}
```

```
##                              label step_major step_minor label_minor
## 5 fit.models_0_Max.cor.Y.Time.Poly          1          4      glmnet
## 6  fit.models_0_Max.cor.Y.Time.Lag          1          5      glmnet
##       bgn     end elapsed
## 5 102.250 114.112  11.863
## 6 114.113      NA      NA
## [1] "fitting model: Max.cor.Y.Time.Lag##rcv#glmnet"
## [1] "    indep_vars: WordCount.root2,NDSS.my.fctr,PubDate.last2.log1p,PubDate.last4.log1p,PubDate.last8.log1p,PubDate.last16.log1p,PubDate.last32.log1p"
## Aggregating results
## Selecting tuning parameters
## Fitting alpha = 0.1, lambda = 0.0934 on full training set
```

```
## Warning in myfit_mdl(mdl_specs_lst = myinit_mdl_specs_lst(mdl_specs_lst
## = list(id.prefix = "Max.cor.Y.Time.Lag", : model's bestTune found at an
## extreme of tuneGrid for parameter: alpha
```

```
## Warning in myfit_mdl(mdl_specs_lst = myinit_mdl_specs_lst(mdl_specs_lst
## = list(id.prefix = "Max.cor.Y.Time.Lag", : model's bestTune found at an
## extreme of tuneGrid for parameter: lambda
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-24.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-25.png) 

```
##             Length Class      Mode     
## a0           100   -none-     numeric  
## beta        2600   dgCMatrix  S4       
## df           100   -none-     numeric  
## dim            2   -none-     numeric  
## lambda       100   -none-     numeric  
## dev.ratio    100   -none-     numeric  
## nulldev        1   -none-     numeric  
## npasses        1   -none-     numeric  
## jerr           1   -none-     numeric  
## offset         1   -none-     logical  
## classnames     2   -none-     character
## call           5   -none-     call     
## nobs           1   -none-     numeric  
## lambdaOpt      1   -none-     numeric  
## xNames        26   -none-     character
## problemType    1   -none-     character
## tuneValue      2   data.frame list     
## obsLevels      2   -none-     character
## [1] "min lambda > lambdaOpt:"
##                         (Intercept)                 NDSS.my.fctr#Mltmd# 
##                        -3.110189977                        -0.122999082 
##           NDSS.my.fctr#Opnn#RmFrDbt        NDSS.my.fctr#Opnn#ThPblcEdtr 
##                        -0.486965112                         1.985537850 
##             NDSS.my.fctr#U.S.#Edctn      NDSS.my.fctrBsnss#BsnssDy#Dlbk 
##                        -0.355905159                        -0.150819398 
## NDSS.my.fctrBsnss#BsnssDy#SmllBsnss     NDSS.my.fctrBsnss#Crsswrds/Gms# 
##                        -0.143861898                         2.409404874 
##              NDSS.my.fctrCltr#Arts#              NDSS.my.fctrFrgn#Wrld# 
##                        -0.135025648                        -0.182763766 
##        NDSS.my.fctrFrgn#Wrld#AsPcfc              NDSS.my.fctrOpEd#Opnn# 
##                        -0.334337141                         2.519450359 
##              NDSS.my.fctrScnc#Hlth#             NDSS.my.fctrStyls##Fshn 
##                         2.011582007                        -0.261047460 
##             NDSS.my.fctrStyls#U.S.#                 NDSS.my.fctrTStyl## 
##                         1.847003574                        -0.424630217 
##              NDSS.my.fctrTrvl#Trvl#                 PubDate.last2.log1p 
##                        -0.108706799                         0.015692430 
##                 PubDate.last4.log1p                 PubDate.last8.log1p 
##                         0.025510119                         0.005128644 
##                     WordCount.root2 
##                         0.031642730 
## [1] "max lambda < lambdaOpt:"
##                         (Intercept)                 NDSS.my.fctr#Mltmd# 
##                        -3.228591381                        -0.157045755 
##           NDSS.my.fctr#Opnn#RmFrDbt        NDSS.my.fctr#Opnn#ThPblcEdtr 
##                        -0.562722304                         2.094599975 
##             NDSS.my.fctr#U.S.#Edctn      NDSS.my.fctrBsnss#BsnssDy#Dlbk 
##                        -0.390420136                        -0.161009275 
## NDSS.my.fctrBsnss#BsnssDy#SmllBsnss     NDSS.my.fctrBsnss#Crsswrds/Gms# 
##                        -0.174871733                         2.498780352 
##              NDSS.my.fctrCltr#Arts#              NDSS.my.fctrFrgn#Wrld# 
##                        -0.143896608                        -0.213094088 
##        NDSS.my.fctrFrgn#Wrld#AsPcfc              NDSS.my.fctrOpEd#Opnn# 
##                        -0.374278240                         2.606761820 
##              NDSS.my.fctrScnc#Hlth#             NDSS.my.fctrStyls##Fshn 
##                         2.089084565                        -0.300271918 
##             NDSS.my.fctrStyls#U.S.#                 NDSS.my.fctrTStyl## 
##                         1.923211176                        -0.451317439 
##              NDSS.my.fctrTrvl#Trvl#                  NDSS.my.fctrmyOthr 
##                        -0.140661020                        -0.021516646 
##                 PubDate.last2.log1p                 PubDate.last4.log1p 
##                         0.017794795                         0.028174208 
##                 PubDate.last8.log1p                     WordCount.root2 
##                         0.007857041                         0.033209968
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-26.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-27.png) 

```
##          Prediction
## Reference    N    Y
##         N 3791  150
##         Y  173  690
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   9.327644e-01   7.694835e-01   9.253120e-01   9.396854e-01   8.203580e-01 
## AccuracyPValue  McnemarPValue 
##  3.123849e-114   2.209097e-01
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-28.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-29.png) 

```
##          Prediction
## Reference   N   Y
##         N 908 590
##         Y  21 209
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   6.464120e-01   2.515036e-01   6.233476e-01   6.689781e-01   8.668981e-01 
## AccuracyPValue  McnemarPValue 
##   1.000000e+00  7.592139e-117 
##                               id
## 1 Max.cor.Y.Time.Lag##rcv#glmnet
##                                                                                                                                feats
## 1 WordCount.root2,NDSS.my.fctr,PubDate.last2.log1p,PubDate.last4.log1p,PubDate.last8.log1p,PubDate.last16.log1p,PubDate.last32.log1p
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1               5                      4.225                 0.318
##   max.AUCpROC.fit max.Sens.fit max.Spec.fit max.AUCROCR.fit
## 1        0.841928    0.9781781    0.7056779       0.9581563
##   opt.prob.threshold.fit max.f.score.fit max.Accuracy.fit
## 1                    0.3       0.8103347        0.9279084
##   max.AccuracyLower.fit max.AccuracyUpper.fit max.Kappa.fit
## 1              0.925312             0.9396854      0.735003
##   max.AUCpROC.OOB max.Sens.OOB max.Spec.OOB max.AUCROCR.OOB
## 1       0.5925379    0.9372497    0.2478261       0.8117635
##   opt.prob.threshold.OOB max.f.score.OOB max.Accuracy.OOB
## 1                    0.1       0.4062196         0.646412
##   max.AccuracyLower.OOB max.AccuracyUpper.OOB max.Kappa.OOB
## 1             0.6233476             0.6689781     0.2515036
##   max.AccuracySD.fit max.KappaSD.fit
## 1        0.004559889      0.02005103
```

```r
# Interactions.High.cor.Y
if (length(int_feats <- setdiff(setdiff(unique(glb_feats_df$cor.high.X), NA), 
                                subset(glb_feats_df, nzv)$id)) > 0) {
    fit.models_0_chunk_df <- myadd_chunk(fit.models_0_chunk_df, 
                    paste0("fit.models_0_", "Interact.High.cor.Y"), major.inc = FALSE,
                                    label.minor = "glmnet")

    ret_lst <- myfit_mdl(mdl_specs_lst=myinit_mdl_specs_lst(mdl_specs_lst=list(
        id.prefix="Interact.High.cor.Y", 
        type=glb_model_type, trainControl.method="repeatedcv",
        trainControl.number=glb_rcv_n_folds, trainControl.repeats=glb_rcv_n_repeats,
            trainControl.classProbs = glb_is_classification,
            trainControl.summaryFunction = glbMdlMetricSummaryFn,
            train.metric = glbMdlMetricSummary, 
            train.maximize = glbMdlMetricMaximize,    
        train.method="glmnet")),
        indep_vars=c(max_cor_y_x_vars, paste(max_cor_y_x_vars[1], int_feats, sep=":")),
        rsp_var=glb_rsp_var, 
        fit_df=glbObsFit, OOB_df=glbObsOOB)
}    
```

```
##                              label step_major step_minor label_minor
## 6  fit.models_0_Max.cor.Y.Time.Lag          1          5      glmnet
## 7 fit.models_0_Interact.High.cor.Y          1          6      glmnet
##       bgn    end elapsed
## 6 114.113 124.99  10.877
## 7 124.990     NA      NA
## [1] "fitting model: Interact.High.cor.Y##rcv#glmnet"
## [1] "    indep_vars: WordCount.root2,NDSS.my.fctr,WordCount.root2:WordCount.root2,WordCount.root2:PubDate.day.minutes.poly.1,WordCount.root2:PubDate.last4.log1p,WordCount.root2:PubDate.month.fctr"
## Aggregating results
## Selecting tuning parameters
## Fitting alpha = 0.775, lambda = 0.000934 on full training set
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-30.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-31.png) 

```
##             Length Class      Mode     
## a0            92   -none-     numeric  
## beta        2392   dgCMatrix  S4       
## df            92   -none-     numeric  
## dim            2   -none-     numeric  
## lambda        92   -none-     numeric  
## dev.ratio     92   -none-     numeric  
## nulldev        1   -none-     numeric  
## npasses        1   -none-     numeric  
## jerr           1   -none-     numeric  
## offset         1   -none-     logical  
## classnames     2   -none-     character
## call           5   -none-     call     
## nobs           1   -none-     numeric  
## lambdaOpt      1   -none-     numeric  
## xNames        26   -none-     character
## problemType    1   -none-     character
## tuneValue      2   data.frame list     
## obsLevels      2   -none-     character
## [1] "min lambda > lambdaOpt:"
##                                (Intercept) 
##                              -4.8684452161 
##                        NDSS.my.fctr#Mltmd# 
##                              -1.0487925341 
##                  NDSS.my.fctr#Opnn#RmFrDbt 
##                              -5.4745572239 
##               NDSS.my.fctr#Opnn#ThPblcEdtr 
##                               4.2885428864 
##                    NDSS.my.fctr#U.S.#Edctn 
##                              -2.8325574053 
##             NDSS.my.fctrBsnss#BsnssDy#Dlbk 
##                              -0.3151450444 
##        NDSS.my.fctrBsnss#BsnssDy#SmllBsnss 
##                              -0.7572259115 
##            NDSS.my.fctrBsnss#Crsswrds/Gms# 
##                               4.3391439837 
##                 NDSS.my.fctrBsnss#Tchnlgy# 
##                               0.8045020529 
##                     NDSS.my.fctrCltr#Arts# 
##                              -0.3229201625 
##                     NDSS.my.fctrFrgn#Wrld# 
##                              -1.8392467451 
##               NDSS.my.fctrFrgn#Wrld#AsPcfc 
##                              -1.8077498032 
##                  NDSS.my.fctrMtr#N.Y./Rgn# 
##                               0.2627263003 
##                     NDSS.my.fctrOpEd#Opnn# 
##                               4.7415088342 
##                     NDSS.my.fctrScnc#Hlth# 
##                               3.7483695286 
##                    NDSS.my.fctrStyls##Fshn 
##                              -2.4194040890 
##                    NDSS.my.fctrStyls#U.S.# 
##                               3.4229733393 
##                        NDSS.my.fctrTStyl## 
##                              -2.0211781889 
##                     NDSS.my.fctrTrvl#Trvl# 
##                              -1.7751305343 
##                         NDSS.my.fctrmyOthr 
##                              -2.3074809466 
##                            WordCount.root2 
##                               0.0361357822 
## WordCount.root2:PubDate.day.minutes.poly.1 
##                               1.0813556024 
##        WordCount.root2:PubDate.last4.log1p 
##                               0.0069832077 
##       WordCount.root2:PubDate.month.fctr10 
##                               0.0039570640 
##       WordCount.root2:PubDate.month.fctr11 
##                              -0.0001256918 
## [1] "max lambda < lambdaOpt:"
##                                (Intercept) 
##                              -4.8709915569 
##                        NDSS.my.fctr#Mltmd# 
##                              -1.0849702043 
##                  NDSS.my.fctr#Opnn#RmFrDbt 
##                              -5.5929746399 
##               NDSS.my.fctr#Opnn#ThPblcEdtr 
##                               4.2928665989 
##                    NDSS.my.fctr#U.S.#Edctn 
##                              -2.9379189459 
##             NDSS.my.fctrBsnss#BsnssDy#Dlbk 
##                              -0.3276698231 
##        NDSS.my.fctrBsnss#BsnssDy#SmllBsnss 
##                              -0.7772032984 
##            NDSS.my.fctrBsnss#Crsswrds/Gms# 
##                               4.3366629069 
##                 NDSS.my.fctrBsnss#Tchnlgy# 
##                               0.8004536454 
##                     NDSS.my.fctrCltr#Arts# 
##                              -0.3372244058 
##                     NDSS.my.fctrFrgn#Wrld# 
##                              -1.9323067428 
##               NDSS.my.fctrFrgn#Wrld#AsPcfc 
##                              -1.8410299268 
##                  NDSS.my.fctrMtr#N.Y./Rgn# 
##                               0.2641549470 
##                     NDSS.my.fctrOpEd#Opnn# 
##                               4.7402435031 
##                     NDSS.my.fctrScnc#Hlth# 
##                               3.7456718341 
##                    NDSS.my.fctrStyls##Fshn 
##                              -2.5160313263 
##                    NDSS.my.fctrStyls#U.S.# 
##                               3.4196115772 
##                        NDSS.my.fctrTStyl## 
##                              -2.0543357159 
##                     NDSS.my.fctrTrvl#Trvl# 
##                              -1.8687591112 
##                         NDSS.my.fctrmyOthr 
##                              -2.4107173389 
##                            WordCount.root2 
##                               0.0361773869 
## WordCount.root2:PubDate.day.minutes.poly.1 
##                               1.0887365081 
##        WordCount.root2:PubDate.last4.log1p 
##                               0.0070269546 
##       WordCount.root2:PubDate.month.fctr10 
##                               0.0039706413 
##       WordCount.root2:PubDate.month.fctr11 
##                              -0.0002339144
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-32.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-33.png) 

```
##          Prediction
## Reference    N    Y
##         N 3787  154
##         Y  173  690
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   9.319317e-01   7.670549e-01   9.244394e-01   9.388938e-01   8.203580e-01 
## AccuracyPValue  McnemarPValue 
##  2.593356e-112   3.195407e-01
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-34.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-35.png) 

```
##          Prediction
## Reference    N    Y
##         N 1164  334
##         Y   71  159
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   7.656250e-01   3.156027e-01   7.449213e-01   7.854227e-01   8.668981e-01 
## AccuracyPValue  McnemarPValue 
##   1.000000e+00   9.555643e-39 
##                                id
## 1 Interact.High.cor.Y##rcv#glmnet
##                                                                                                                                                                            feats
## 1 WordCount.root2,NDSS.my.fctr,WordCount.root2:WordCount.root2,WordCount.root2:PubDate.day.minutes.poly.1,WordCount.root2:PubDate.last4.log1p,WordCount.root2:PubDate.month.fctr
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1              25                      5.157                 0.327
##   max.AUCpROC.fit max.Sens.fit max.Spec.fit max.AUCROCR.fit
## 1       0.8776419    0.9626998     0.792584       0.9625372
##   opt.prob.threshold.fit max.f.score.fit max.Accuracy.fit
## 1                    0.4       0.8084359         0.931585
##   max.AccuracyLower.fit max.AccuracyUpper.fit max.Kappa.fit
## 1             0.9244394             0.9388938      0.764104
##   max.AUCpROC.OOB max.Sens.OOB max.Spec.OOB max.AUCROCR.OOB
## 1       0.6009259    0.9105474    0.2913043       0.8140971
##   opt.prob.threshold.OOB max.f.score.OOB max.Accuracy.OOB
## 1                    0.1        0.439834         0.765625
##   max.AccuracyLower.OOB max.AccuracyUpper.OOB max.Kappa.OOB
## 1             0.7449213             0.7854227     0.3156027
##   max.AccuracySD.fit max.KappaSD.fit
## 1        0.005250654      0.01810996
```

```r
# Low.cor.X
fit.models_0_chunk_df <- myadd_chunk(fit.models_0_chunk_df, 
                        paste0("fit.models_0_", "Low.cor.X"), major.inc = FALSE,
                                     label.minor = "glmnet")
```

```
##                              label step_major step_minor label_minor
## 7 fit.models_0_Interact.High.cor.Y          1          6      glmnet
## 8           fit.models_0_Low.cor.X          1          7      glmnet
##       bgn     end elapsed
## 7 124.990 136.861  11.871
## 8 136.861      NA      NA
```

```r
indep_vars <- subset(glb_feats_df, is.na(cor.high.X) & !nzv & 
                              (exclude.as.feat != 1))[, "id"]  
indep_vars <- myadjust_interaction_feats(indep_vars)
ret_lst <- myfit_mdl(mdl_specs_lst=myinit_mdl_specs_lst(mdl_specs_lst=list(
        id.prefix="Low.cor.X", 
        type=glb_model_type, 
        tune.df = glbMdlTuneParams,        
        trainControl.method="repeatedcv",
        trainControl.number=glb_rcv_n_folds, trainControl.repeats=glb_rcv_n_repeats,
            trainControl.classProbs = glb_is_classification,
            trainControl.summaryFunction = glbMdlMetricSummaryFn,
            train.metric = glbMdlMetricSummary, 
            train.maximize = glbMdlMetricMaximize,    
        train.method="glmnet")),
        indep_vars=indep_vars, rsp_var=glb_rsp_var, 
        fit_df=glbObsFit, OOB_df=glbObsOOB)
```

```
## [1] "fitting model: Low.cor.X##rcv#glmnet"
## [1] "    indep_vars: WordCount.root2,NDSS.my.fctr,PubDate.day.minutes.poly.1,PubDate.wkend,PubDate.day.minutes.poly.4,PubDate.day.minutes.poly.2,PubDate.last4.log1p,PubDate.last16.log1p,PubDate.day.minutes.poly.3,PubDate.month.fctr,.rnorm,PubDate.last32.log1p,PubDate.date.fctr,PubDate.second.fctr,PubDate.minute.fctr,PubDate.wkday.fctr,WordCount.nexp,PubDate.day.minutes.poly.5"
## Aggregating results
## Selecting tuning parameters
## Fitting alpha = 0.1, lambda = 0.0201 on full training set
```

```
## Warning in myfit_mdl(mdl_specs_lst = myinit_mdl_specs_lst(mdl_specs_lst
## = list(id.prefix = "Low.cor.X", : model's bestTune found at an extreme of
## tuneGrid for parameter: alpha
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-36.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-37.png) 

```
##             Length Class      Mode     
## a0           100   -none-     numeric  
## beta        5100   dgCMatrix  S4       
## df           100   -none-     numeric  
## dim            2   -none-     numeric  
## lambda       100   -none-     numeric  
## dev.ratio    100   -none-     numeric  
## nulldev        1   -none-     numeric  
## npasses        1   -none-     numeric  
## jerr           1   -none-     numeric  
## offset         1   -none-     logical  
## classnames     2   -none-     character
## call           5   -none-     call     
## nobs           1   -none-     numeric  
## lambdaOpt      1   -none-     numeric  
## xNames        51   -none-     character
## problemType    1   -none-     character
## tuneValue      2   data.frame list     
## obsLevels      2   -none-     character
## [1] "min lambda > lambdaOpt:"
##                         (Intercept)                 NDSS.my.fctr#Mltmd# 
##                         -4.62776335                         -0.55058989 
##           NDSS.my.fctr#Opnn#RmFrDbt        NDSS.my.fctr#Opnn#ThPblcEdtr 
##                         -2.17201828                          3.43746156 
##             NDSS.my.fctr#U.S.#Edctn      NDSS.my.fctrBsnss#BsnssDy#Dlbk 
##                         -0.90024179                         -0.25001208 
## NDSS.my.fctrBsnss#BsnssDy#SmllBsnss     NDSS.my.fctrBsnss#Crsswrds/Gms# 
##                         -0.51959279                          3.40386794 
##          NDSS.my.fctrBsnss#Tchnlgy#              NDSS.my.fctrCltr#Arts# 
##                          0.41727883                         -0.35499218 
##              NDSS.my.fctrFrgn#Wrld#        NDSS.my.fctrFrgn#Wrld#AsPcfc 
##                         -0.71609323                         -0.95587658 
##              NDSS.my.fctrOpEd#Opnn#              NDSS.my.fctrScnc#Hlth# 
##                          3.79502666                          3.05301675 
##             NDSS.my.fctrStyls##Fshn             NDSS.my.fctrStyls#U.S.# 
##                         -0.97806123                          2.87280463 
##                 NDSS.my.fctrTStyl##              NDSS.my.fctrTrvl#Trvl# 
##                         -0.92999674                         -0.72203823 
##                  NDSS.my.fctrmyOthr             PubDate.date.fctr(7,13] 
##                         -0.87859009                          0.05462536 
##            PubDate.date.fctr(13,19]            PubDate.date.fctr(25,31] 
##                         -0.04859970                          0.03936760 
##          PubDate.day.minutes.poly.1          PubDate.day.minutes.poly.2 
##                         16.52129547                          8.15155112 
##          PubDate.day.minutes.poly.4                PubDate.last16.log1p 
##                          7.63942375                          0.03904314 
##                 PubDate.last4.log1p      PubDate.minute.fctr(14.8,29.5] 
##                          0.05261649                          0.04081341 
##      PubDate.minute.fctr(29.5,44.2]      PubDate.minute.fctr(44.2,59.1] 
##                         -0.06247733                          0.15574933 
##                PubDate.month.fctr10                PubDate.month.fctr11 
##                          0.01965406                         -0.04371938 
##      PubDate.second.fctr(29.5,44.2]                 PubDate.wkday.fctr1 
##                          0.02855656                          0.15316236 
##                 PubDate.wkday.fctr2                 PubDate.wkday.fctr4 
##                         -0.04811783                         -0.02920222 
##                 PubDate.wkday.fctr5                 PubDate.wkday.fctr6 
##                         -0.11794295                         -0.06053284 
##                       PubDate.wkend                      WordCount.nexp 
##                          0.30065231                         -0.01458288 
##                     WordCount.root2 
##                          0.05980150 
## [1] "max lambda < lambdaOpt:"
##                         (Intercept)                 NDSS.my.fctr#Mltmd# 
##                        -4.735626979                        -0.580286171 
##           NDSS.my.fctr#Opnn#RmFrDbt        NDSS.my.fctr#Opnn#ThPblcEdtr 
##                        -2.278694860                         3.498981284 
##             NDSS.my.fctr#U.S.#Edctn      NDSS.my.fctrBsnss#BsnssDy#Dlbk 
##                        -0.945806825                        -0.253070361 
## NDSS.my.fctrBsnss#BsnssDy#SmllBsnss     NDSS.my.fctrBsnss#Crsswrds/Gms# 
##                        -0.538161808                         3.460877704 
##          NDSS.my.fctrBsnss#Tchnlgy#              NDSS.my.fctrCltr#Arts# 
##                         0.443459637                        -0.362873872 
##              NDSS.my.fctrFrgn#Wrld#        NDSS.my.fctrFrgn#Wrld#AsPcfc 
##                        -0.755050599                        -0.999081408 
##              NDSS.my.fctrOpEd#Opnn#              NDSS.my.fctrScnc#Hlth# 
##                         3.856808591                         3.096140054 
##             NDSS.my.fctrStyls##Fshn             NDSS.my.fctrStyls#U.S.# 
##                        -1.029289668                         2.914022368 
##                 NDSS.my.fctrTStyl##              NDSS.my.fctrTrvl#Trvl# 
##                        -0.962511600                        -0.760297888 
##                  NDSS.my.fctrmyOthr             PubDate.date.fctr(7,13] 
##                        -0.939702864                         0.060340544 
##            PubDate.date.fctr(13,19]            PubDate.date.fctr(25,31] 
##                        -0.053309854                         0.046622969 
##          PubDate.day.minutes.poly.1          PubDate.day.minutes.poly.2 
##                        16.771123786                         8.372523804 
##          PubDate.day.minutes.poly.4                PubDate.last16.log1p 
##                         7.771940232                         0.043390584 
##                 PubDate.last4.log1p      PubDate.minute.fctr(14.8,29.5] 
##                         0.053331182                         0.051567949 
##      PubDate.minute.fctr(29.5,44.2]      PubDate.minute.fctr(44.2,59.1] 
##                        -0.063196836                         0.170297377 
##                PubDate.month.fctr10                PubDate.month.fctr11 
##                         0.021788115                        -0.051764776 
##      PubDate.second.fctr(29.5,44.2]      PubDate.second.fctr(44.2,59.1] 
##                         0.035204623                        -0.002666142 
##                 PubDate.wkday.fctr1                 PubDate.wkday.fctr2 
##                         0.160804999                        -0.055611629 
##                 PubDate.wkday.fctr4                 PubDate.wkday.fctr5 
##                        -0.036954495                        -0.127514309 
##                 PubDate.wkday.fctr6                       PubDate.wkend 
##                        -0.089576187                         0.309057681 
##                      WordCount.nexp                     WordCount.root2 
##                        -0.032095197                         0.061485177
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-38.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-39.png) 

```
##          Prediction
## Reference    N    Y
##         N 3743  198
##         Y  138  725
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   9.300583e-01   7.689739e-01   9.224771e-01   9.371115e-01   8.203580e-01 
## AccuracyPValue  McnemarPValue 
##  4.458834e-108   1.287669e-03
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-40.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-41.png) 

```
##          Prediction
## Reference    N    Y
##         N 1136  362
##         Y   60  170
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   7.557870e-01   3.197713e-01   7.348172e-01   7.758846e-01   8.668981e-01 
## AccuracyPValue  McnemarPValue 
##   1.000000e+00   1.299428e-48 
##                      id
## 1 Low.cor.X##rcv#glmnet
##                                                                                                                                                                                                                                                                                                                                                                   feats
## 1 WordCount.root2,NDSS.my.fctr,PubDate.day.minutes.poly.1,PubDate.wkend,PubDate.day.minutes.poly.4,PubDate.day.minutes.poly.2,PubDate.last4.log1p,PubDate.last16.log1p,PubDate.day.minutes.poly.3,PubDate.month.fctr,.rnorm,PubDate.last32.log1p,PubDate.date.fctr,PubDate.second.fctr,PubDate.minute.fctr,PubDate.wkday.fctr,WordCount.nexp,PubDate.day.minutes.poly.5
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1              25                      6.823                 0.515
##   max.AUCpROC.fit max.Sens.fit max.Spec.fit max.AUCROCR.fit
## 1       0.8751806    0.9647298    0.7856315       0.9622982
##   opt.prob.threshold.fit max.f.score.fit max.Accuracy.fit
## 1                    0.2       0.8118701        0.9323486
##   max.AccuracyLower.fit max.AccuracyUpper.fit max.Kappa.fit
## 1             0.9224771             0.9371115     0.7651127
##   max.AUCpROC.OOB max.Sens.OOB max.Spec.OOB max.AUCROCR.OOB
## 1       0.6015934    0.9118825    0.2913043       0.8160359
##   opt.prob.threshold.OOB max.f.score.OOB max.Accuracy.OOB
## 1                    0.1       0.4461942         0.755787
##   max.AccuracyLower.OOB max.AccuracyUpper.OOB max.Kappa.OOB
## 1             0.7348172             0.7758846     0.3197713
##   max.AccuracySD.fit max.KappaSD.fit
## 1        0.005383249       0.0186909
```

```r
fit.models_0_chunk_df <- 
    myadd_chunk(fit.models_0_chunk_df, "fit.models_0_end", major.inc = FALSE,
                label.minor = "teardown")
```

```
##                    label step_major step_minor label_minor     bgn     end
## 8 fit.models_0_Low.cor.X          1          7      glmnet 136.861 150.646
## 9       fit.models_0_end          1          8    teardown 150.646      NA
##   elapsed
## 8  13.785
## 9      NA
```

```r
rm(ret_lst)

glb_chunks_df <- myadd_chunk(glb_chunks_df, "fit.models", major.inc = FALSE)
```

```
##         label step_major step_minor label_minor    bgn     end elapsed
## 10 fit.models          6          0           0  80.36 150.659  70.299
## 11 fit.models          6          1           1 150.66      NA      NA
```


```r
fit.models_1_chunk_df <- myadd_chunk(NULL, "fit.models_1_bgn", label.minor="setup")
```

```
##              label step_major step_minor label_minor     bgn end elapsed
## 1 fit.models_1_bgn          1          0       setup 156.765  NA      NA
```

```r
#stop(here"); glb2Sav(); all.equal(glb_models_df, sav_models_df)
topindep_var <- NULL; interact_vars <- NULL;
for (mdl_id_pfx in names(glb_mdl_family_lst)) {
    fit.models_1_chunk_df <- 
        myadd_chunk(fit.models_1_chunk_df, paste0("fit.models_1_", mdl_id_pfx),
                    major.inc = FALSE, label.minor = "setup")

    indep_vars <- NULL;

    if (grepl("\\.Interact", mdl_id_pfx)) {
        if (is.null(topindep_var) && is.null(interact_vars)) {
        #   select best glmnet model upto now
            dsp_models_df <- orderBy(model_sel_frmla <- get_model_sel_frmla(),
                                     glb_models_df)
            dsp_models_df <- subset(dsp_models_df, 
                                    grepl(".glmnet", id, fixed = TRUE))
            bst_mdl_id <- dsp_models_df$id[1]
            mdl_id_pfx <- 
                paste(c(head(unlist(strsplit(bst_mdl_id, "[.]")), -1), "Interact"),
                      collapse=".")
        #   select important features
            if (is.null(bst_featsimp_df <- 
                        myget_feats_importance(glb_models_lst[[bst_mdl_id]]))) {
                warning("Base model for RFE.Interact: ", bst_mdl_id, 
                        " has no important features")
                next
            }    
            
            topindep_ix <- 1
            while (is.null(topindep_var) && (topindep_ix <= nrow(bst_featsimp_df))) {
                topindep_var <- row.names(bst_featsimp_df)[topindep_ix]
                if (grepl(".fctr", topindep_var, fixed=TRUE))
                    topindep_var <- 
                        paste0(unlist(strsplit(topindep_var, ".fctr"))[1], ".fctr")
                if (topindep_var %in% names(glbFeatsInteractionOnly)) {
                    topindep_var <- NULL; topindep_ix <- topindep_ix + 1
                } else break
            }
            
        #   select features with importance > max(10, importance of .rnorm) & is not highest
        #       combine factor dummy features to just the factor feature
            if (length(pos_rnorm <- 
                       grep(".rnorm", row.names(bst_featsimp_df), fixed=TRUE)) > 0)
                imp_rnorm <- bst_featsimp_df[pos_rnorm, 1] else
                imp_rnorm <- NA    
            imp_cutoff <- max(10, imp_rnorm, na.rm=TRUE)
            interact_vars <- 
                tail(row.names(subset(bst_featsimp_df, 
                                      imp > imp_cutoff)), -1)
            if (length(interact_vars) > 0) {
                interact_vars <-
                    myadjust_interaction_feats(myextract_actual_feats(interact_vars))
                interact_vars <- 
                    interact_vars[!grepl(topindep_var, interact_vars, fixed=TRUE)]
            }
            ### bid0_sp only
#             interact_vars <- c(
#     "biddable", "D.ratio.sum.TfIdf.wrds.n", "D.TfIdf.sum.stem.stop.Ratio", "D.sum.TfIdf",
#     "D.TfIdf.sum.post.stop", "D.TfIdf.sum.post.stem", "D.ratio.wrds.stop.n.wrds.n", "D.chrs.uppr.n.log",
#     "D.chrs.n.log", "color.fctr"
#     # , "condition.fctr", "prdl.my.descr.fctr"
#                                 )
#            interact_vars <- setdiff(interact_vars, c("startprice.dgt2.is9", "color.fctr"))
            ###
            indep_vars <- myextract_actual_feats(row.names(bst_featsimp_df))
            indep_vars <- setdiff(indep_vars, topindep_var)
            if (length(interact_vars) > 0) {
                indep_vars <- 
                    setdiff(indep_vars, myextract_actual_feats(interact_vars))
                indep_vars <- c(indep_vars, 
                    paste(topindep_var, setdiff(interact_vars, topindep_var), 
                          sep = "*"))
            } else indep_vars <- union(indep_vars, topindep_var)
        }
    }
    
    if (is.null(indep_vars))
        indep_vars <- glb_mdl_feats_lst[[mdl_id_pfx]]

    if (is.null(indep_vars) && grepl("RFE\\.", mdl_id_pfx))
        indep_vars <- myextract_actual_feats(predictors(rfe_fit_results))
    
    if (is.null(indep_vars))
        indep_vars <- subset(glb_feats_df, !nzv & (exclude.as.feat != 1))[, "id"]
    
    if ((length(indep_vars) == 1) && (grepl("^%<d-%", indep_vars))) {    
        indep_vars <- 
            eval(parse(text = str_trim(unlist(strsplit(indep_vars, "%<d-%"))[2])))
    }    

    indep_vars <- myadjust_interaction_feats(indep_vars)
    
    if (grepl("\\.Interact", mdl_id_pfx)) { 
        # if (method != tail(unlist(strsplit(bst_mdl_id, "[.]")), 1)) next
        if (is.null(glb_mdl_family_lst[[mdl_id_pfx]])) {
            if (!is.null(glb_mdl_family_lst[["Best.Interact"]]))
                glb_mdl_family_lst[[mdl_id_pfx]] <-
                    glb_mdl_family_lst[["Best.Interact"]]
        }
    }
    
    if (!is.null(glbObsFitOutliers[[mdl_id_pfx]])) {
        fitobs_df <- glbObsFit[!(glbObsFit[, glb_id_var] %in%
                                         glbObsFitOutliers[[mdl_id_pfx]]), ]
    } else fitobs_df <- glbObsFit

    if (is.null(glb_mdl_family_lst[[mdl_id_pfx]]))
        mdl_methods <- glbMdlMethods else
        mdl_methods <- glb_mdl_family_lst[[mdl_id_pfx]]    

    for (method in mdl_methods) {
        if (method %in% c("rpart", "rf")) {
            # rpart:    fubar's the tree
            # rf:       skip the scenario w/ .rnorm for speed
            indep_vars <- setdiff(indep_vars, c(".rnorm"))
            #mdl_id <- paste0(mdl_id_pfx, ".no.rnorm")
        } 

        fit.models_1_chunk_df <- myadd_chunk(fit.models_1_chunk_df, 
                            paste0("fit.models_1_", mdl_id_pfx), major.inc = FALSE,
                                    label.minor = method)
# The last([[:digit:]]+)(.*)\\.ctg feats are taking a long time for this experiment & nor proving to be important
indep_vars <- indep_vars[!grepl("\\.last([[:digit:]]+)(.*)\\.ctg", indep_vars)]
# The poly\\.([[:digit:]]+)\\.ctg feats are taking a long time for this experiment & nor proving to be important
indep_vars <- indep_vars[!grepl("\\.poly\\.([[:digit:]]+)\\.ctg", indep_vars)]

        ret_lst <- 
            myfit_mdl(mdl_specs_lst = myinit_mdl_specs_lst(mdl_specs_lst = list(
            id.prefix = mdl_id_pfx, 
            type = glb_model_type, 
            tune.df = glbMdlTuneParams,
            trainControl.method = "repeatedcv",
            trainControl.number = glb_rcv_n_folds,
            trainControl.repeats = glb_rcv_n_repeats,
            trainControl.classProbs = glb_is_classification,
            trainControl.summaryFunction = glbMdlMetricSummaryFn,
#trainControl.allowParallel = FALSE,
            train.metric = glbMdlMetricSummary, 
            train.maximize = glbMdlMetricMaximize,    
            train.method = method)),
            indep_vars = indep_vars, rsp_var = glb_rsp_var, 
            fit_df = fitobs_df, OOB_df = glbObsOOB)
    }
}
```

```
##                label step_major step_minor label_minor     bgn     end
## 1   fit.models_1_bgn          1          0       setup 156.765 156.775
## 2 fit.models_1_All.X          1          1       setup 156.775      NA
##   elapsed
## 1    0.01
## 2      NA
##                label step_major step_minor label_minor     bgn     end
## 2 fit.models_1_All.X          1          1       setup 156.775 156.781
## 3 fit.models_1_All.X          1          2      glmnet 156.782      NA
##   elapsed
## 2   0.006
## 3      NA
## [1] "fitting model: All.X##rcv#glmnet"
## [1] "    indep_vars: WordCount.root2,WordCount.log1p,NDSS.my.fctr,PubDate.day.minutes.poly.1,PubDate.hour.fctr,PubDate.wkend,PubDate.day.minutes.poly.4,PubDate.day.minutes.poly.2,PubDate.last4.log1p,PubDate.last2.log1p,PubDate.last8.log1p,PubDate.last16.log1p,PubDate.day.minutes.poly.3,PubDate.month.fctr,PubDate.juliandate,.rnorm,PubDate.last32.log1p,PubDate.date.fctr,PubDate.second.fctr,PubDate.minute.fctr,PubDate.wkday.fctr,WordCount.nexp,PubDate.day.minutes.poly.5"
## Aggregating results
## Selecting tuning parameters
## Fitting alpha = 0.775, lambda = 0.0201 on full training set
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_1-1.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_1-2.png) 

```
##             Length Class      Mode     
## a0            89   -none-     numeric  
## beta        5073   dgCMatrix  S4       
## df            89   -none-     numeric  
## dim            2   -none-     numeric  
## lambda        89   -none-     numeric  
## dev.ratio     89   -none-     numeric  
## nulldev        1   -none-     numeric  
## npasses        1   -none-     numeric  
## jerr           1   -none-     numeric  
## offset         1   -none-     logical  
## classnames     2   -none-     character
## call           5   -none-     call     
## nobs           1   -none-     numeric  
## lambdaOpt      1   -none-     numeric  
## xNames        57   -none-     character
## problemType    1   -none-     character
## tuneValue      2   data.frame list     
## obsLevels      2   -none-     character
## [1] "min lambda > lambdaOpt:"
##                     (Intercept)       NDSS.my.fctr#Opnn#RmFrDbt 
##                     -4.15778549                     -0.29068581 
##    NDSS.my.fctr#Opnn#ThPblcEdtr NDSS.my.fctrBsnss#Crsswrds/Gms# 
##                      2.70584441                      3.68584256 
##          NDSS.my.fctrOpEd#Opnn#          NDSS.my.fctrScnc#Hlth# 
##                      4.00374690                      3.09723359 
##         NDSS.my.fctrStyls#U.S.#             NDSS.my.fctrTStyl## 
##                      2.85084048                     -0.13236677 
##      PubDate.day.minutes.poly.1                 WordCount.log1p 
##                      5.59078306                      0.08550402 
##                 WordCount.root2 
##                      0.04496835 
## [1] "max lambda < lambdaOpt:"
##                     (Intercept)       NDSS.my.fctr#Opnn#RmFrDbt 
##                    -4.284190646                    -0.447707141 
##    NDSS.my.fctr#Opnn#ThPblcEdtr NDSS.my.fctrBsnss#Crsswrds/Gms# 
##                     2.827765103                     3.751516261 
##      NDSS.my.fctrBsnss#Tchnlgy#          NDSS.my.fctrOpEd#Opnn# 
##                     0.081535690                     4.070511321 
##          NDSS.my.fctrScnc#Hlth#         NDSS.my.fctrStyls#U.S.# 
##                     3.163208431                     2.914061199 
##             NDSS.my.fctrTStyl##      PubDate.day.minutes.poly.1 
##                    -0.187607543                     7.038308760 
##                   PubDate.wkend                 WordCount.log1p 
##                     0.008146344                     0.096036181 
##                 WordCount.root2 
##                     0.046479071
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_1-3.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_1-4.png) 

```
##          Prediction
## Reference    N    Y
##         N 3798  143
##         Y  178  685
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   9.331807e-01   7.696469e-01   9.257484e-01   9.400811e-01   8.203580e-01 
## AccuracyPValue  McnemarPValue 
##  3.361187e-115   5.773628e-02
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_1-5.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_1-6.png) 

```
##          Prediction
## Reference    N    Y
##         N 1195  303
##         Y   87  143
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   7.743056e-01   3.001637e-01   7.538491e-01   7.938262e-01   8.668981e-01 
## AccuracyPValue  McnemarPValue 
##   1.000000e+00   1.330234e-27 
##                  id
## 1 All.X##rcv#glmnet
##                                                                                                                                                                                                                                                                                                                                                                                                                                                                feats
## 1 WordCount.root2,WordCount.log1p,NDSS.my.fctr,PubDate.day.minutes.poly.1,PubDate.hour.fctr,PubDate.wkend,PubDate.day.minutes.poly.4,PubDate.day.minutes.poly.2,PubDate.last4.log1p,PubDate.last2.log1p,PubDate.last8.log1p,PubDate.last16.log1p,PubDate.day.minutes.poly.3,PubDate.month.fctr,PubDate.juliandate,.rnorm,PubDate.last32.log1p,PubDate.date.fctr,PubDate.second.fctr,PubDate.minute.fctr,PubDate.wkday.fctr,WordCount.nexp,PubDate.day.minutes.poly.5
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1              25                      9.521                 0.671
##   max.AUCpROC.fit max.Sens.fit max.Spec.fit max.AUCROCR.fit
## 1       0.8734975    0.9659985    0.7809965       0.9528396
##   opt.prob.threshold.fit max.f.score.fit max.Accuracy.fit
## 1                    0.4       0.8101715        0.9320707
##   max.AccuracyLower.fit max.AccuracyUpper.fit max.Kappa.fit
## 1             0.9257484             0.9400811      0.763027
##   max.AUCpROC.OOB max.Sens.OOB max.Spec.OOB max.AUCROCR.OOB
## 1       0.5950717    0.9118825    0.2782609       0.7995951
##   opt.prob.threshold.OOB max.f.score.OOB max.Accuracy.OOB
## 1                    0.1       0.4230769        0.7743056
##   max.AccuracyLower.OOB max.AccuracyUpper.OOB max.Kappa.OOB
## 1             0.7538491             0.7938262     0.3001637
##   max.AccuracySD.fit max.KappaSD.fit
## 1        0.005666511      0.02009517
```

```r
# Check if other preProcess methods improve model performance
fit.models_1_chunk_df <- 
    myadd_chunk(fit.models_1_chunk_df, "fit.models_1_preProc", major.inc = FALSE,
                label.minor = "preProc")
```

```
##                  label step_major step_minor label_minor     bgn     end
## 3   fit.models_1_All.X          1          2      glmnet 156.782 174.144
## 4 fit.models_1_preProc          1          3     preProc 174.144      NA
##   elapsed
## 3  17.362
## 4      NA
```

```r
mdl_id <- orderBy(get_model_sel_frmla(), glb_models_df)[1, "id"]
indep_vars_vctr <- trim(unlist(strsplit(glb_models_df[glb_models_df$id == mdl_id,
                                                      "feats"], "[,]")))
method <- tail(unlist(strsplit(mdl_id, "[.]")), 1)
mdl_id_pfx <- paste0(head(unlist(strsplit(mdl_id, "[.]")), -1), collapse = ".")
if (!is.null(glbObsFitOutliers[[mdl_id_pfx]])) {
    fitobs_df <- glbObsFit[!(glbObsFit[, glb_id_var] %in%
                                     glbObsFitOutliers[[mdl_id_pfx]]), ]
} else fitobs_df <- glbObsFit

for (prePr in glb_preproc_methods) {   
    # The operations are applied in this order: 
    #   Box-Cox/Yeo-Johnson transformation, centering, scaling, range, imputation, PCA, ICA then spatial sign.
    
    ret_lst <- myfit_mdl(mdl_specs_lst=myinit_mdl_specs_lst(mdl_specs_lst=list(
            id.prefix=mdl_id_pfx, 
            type=glb_model_type, tune.df=glbMdlTuneParams,
            trainControl.method="repeatedcv",
            trainControl.number=glb_rcv_n_folds,
            trainControl.repeats=glb_rcv_n_repeats,
            trainControl.classProbs = glb_is_classification,
            trainControl.summaryFunction = glbMdlMetricSummaryFn,
            train.metric = glbMdlMetricSummary, 
            train.maximize = glbMdlMetricMaximize,    
            train.method=method, train.preProcess=prePr)),
            indep_vars=indep_vars_vctr, rsp_var=glb_rsp_var, 
            fit_df=fitobs_df, OOB_df=glbObsOOB)
}            
    
    # If (All|RFE).X.glm is less accurate than Low.Cor.X.glm
    #   check NA coefficients & filter appropriate terms in indep_vars_vctr
#     if (method == "glm") {
#         orig_glm <- glb_models_lst[[paste0(mdl_id, ".", model_method)]]$finalModel
#         orig_glm <- glb_models_lst[["All.X.glm"]]$finalModel; print(summary(orig_glm))
#         orig_glm <- glb_models_lst[["RFE.X.glm"]]$finalModel; print(summary(orig_glm))
#           require(car)
#           vif_orig_glm <- vif(orig_glm); print(vif_orig_glm)
#           # if vif errors out with "there are aliased coefficients in the model"
#               alias_orig_glm <- alias(orig_glm); alias_complete_orig_glm <- (alias_orig_glm$Complete > 0); alias_complete_orig_glm <- alias_complete_orig_glm[rowSums(alias_complete_orig_glm) > 0, colSums(alias_complete_orig_glm) > 0]; print(alias_complete_orig_glm)
#           print(vif_orig_glm[!is.na(vif_orig_glm) & (vif_orig_glm == Inf)])
#           print(which.max(vif_orig_glm))
#           print(sort(vif_orig_glm[vif_orig_glm >= 1.0e+03], decreasing=TRUE))
#           glbObsFit[c(1143, 3637, 3953, 4105), c("UniqueID", "Popular", "H.P.quandary", "Headline")]
#           glb_feats_df[glb_feats_df$id %in% grep("[HSA]\\.chrs.n.log", glb_feats_df$id, value=TRUE) | glb_feats_df$cor.high.X %in%    grep("[HSA]\\.chrs.n.log", glb_feats_df$id, value=TRUE), ]
#           all.equal(glbObsAll$S.chrs.uppr.n.log, glbObsAll$A.chrs.uppr.n.log)
#           cor(glbObsAll$S.T.herald, glbObsAll$S.T.tribun)
#           mydspObs(Abstract.contains="[Dd]iar", cols=("Abstract"), all=TRUE)
#           subset(glb_feats_df, cor.y.abs <= glb_feats_df[glb_feats_df$id == ".rnorm", "cor.y.abs"])
#         corxx_mtrx <- cor(data.matrix(glbObsAll[, setdiff(names(glbObsAll), myfind_chr_cols_df(glbObsAll))]), use="pairwise.complete.obs"); abs_corxx_mtrx <- abs(corxx_mtrx); diag(abs_corxx_mtrx) <- 0
#           which.max(abs_corxx_mtrx["S.T.tribun", ])
#           abs_corxx_mtrx["A.npnct08.log", "S.npnct08.log"]
#         step_glm <- step(orig_glm)
#     }
    # Since caret does not optimize rpart well
#     if (method == "rpart")
#         ret_lst <- myfit_mdl(mdl_id=paste0(mdl_id_pfx, ".cp.0"), model_method=method,
#                                 indep_vars_vctr=indep_vars_vctr,
#                                 model_type=glb_model_type,
#                                 rsp_var=glb_rsp_var,
#                                 fit_df=glbObsFit, OOB_df=glbObsOOB,        
#             n_cv_folds=0, tune_models_df=data.frame(parameter="cp", min=0.0, max=0.0, by=0.1))

# User specified
#   Ensure at least 2 vars in each regression; else varImp crashes
# sav_models_lst <- glb_models_lst; sav_models_df <- glb_models_df; sav_featsimp_df <- glb_featsimp_df; all.equal(sav_featsimp_df, glb_featsimp_df)
# glb_models_lst <- sav_models_lst; glb_models_df <- sav_models_df; glm_featsimp_df <- sav_featsimp_df

    # easier to exclude features
# require(gdata) # needed for trim
# mdl_id <- "";
# indep_vars_vctr <- head(subset(glb_models_df, grepl("All\\.X\\.", mdl_id), select=feats)
#                         , 1)[, "feats"]
# indep_vars_vctr <- trim(unlist(strsplit(indep_vars_vctr, "[,]")))
# indep_vars_vctr <- setdiff(indep_vars_vctr, ".rnorm")

    # easier to include features
#stop(here"); sav_models_df <- glb_models_df; glb_models_df <- sav_models_df
# !_sp
# mdl_id <- "csm"; indep_vars_vctr <- c(NULL
#     ,"prdline.my.fctr", "prdline.my.fctr:.clusterid.fctr"
#     ,"prdline.my.fctr*biddable"
#     #,"prdline.my.fctr*startprice.log"
#     #,"prdline.my.fctr*startprice.diff"    
#     ,"prdline.my.fctr*condition.fctr"
#     ,"prdline.my.fctr*D.terms.post.stop.n"
#     #,"prdline.my.fctr*D.terms.post.stem.n"
#     ,"prdline.my.fctr*cellular.fctr"    
# #    ,"<feat1>:<feat2>"
#                                            )
# for (method in glbMdlMethods) {
#     ret_lst <- myfit_mdl(mdl_id=mdl_id, model_method=method,
#                                 indep_vars_vctr=indep_vars_vctr,
#                                 model_type=glb_model_type,
#                                 rsp_var=glb_rsp_var,
#                                 fit_df=glbObsFit, OOB_df=glbObsOOB,
#                     n_cv_folds=glb_rcv_n_folds, tune_models_df=glbMdlTuneParams)
#     csm_mdl_id <- paste0(mdl_id, ".", method)
#     csm_featsimp_df <- myget_feats_importance(glb_models_lst[[paste0(mdl_id, ".",
#                                                                      method)]]);               print(head(csm_featsimp_df))
# }
###

# Ntv.1.lm <- lm(reformulate(indep_vars_vctr, glb_rsp_var), glbObsTrn); print(summary(Ntv.1.lm))

#glb_models_df[, "max.Accuracy.OOB", FALSE]
#varImp(glb_models_lst[["Low.cor.X.glm"]])
#orderBy(~ -Overall, varImp(glb_models_lst[["All.X.2.glm"]])$imp)
#orderBy(~ -Overall, varImp(glb_models_lst[["All.X.3.glm"]])$imp)
#glb_feats_df[grepl("npnct28", glb_feats_df$id), ]

    # User specified bivariate models
#     indep_vars_vctr_lst <- list()
#     for (feat in setdiff(names(glbObsFit), 
#                          union(glb_rsp_var, glbFeatsExclude)))
#         indep_vars_vctr_lst[["feat"]] <- feat

    # User specified combinatorial models
#     indep_vars_vctr_lst <- list()
#     combn_mtrx <- combn(c("<feat1_name>", "<feat2_name>", "<featn_name>"), 
#                           <num_feats_to_choose>)
#     for (combn_ix in 1:ncol(combn_mtrx))
#         #print(combn_mtrx[, combn_ix])
#         indep_vars_vctr_lst[[combn_ix]] <- combn_mtrx[, combn_ix]
    
    # template for myfit_mdl
    #   rf is hard-coded in caret to recognize only Accuracy / Kappa evaluation metrics
    #       only for OOB in trainControl ?
    
#     ret_lst <- myfit_mdl_fn(mdl_id=paste0(mdl_id_pfx, ""), model_method=method,
#                             indep_vars_vctr=indep_vars_vctr,
#                             rsp_var=glb_rsp_var,
#                             fit_df=glbObsFit, OOB_df=glbObsOOB,
#                             n_cv_folds=glb_rcv_n_folds, tune_models_df=glbMdlTuneParams,
#                             model_loss_mtrx=glbMdlMetric_terms,
#                             model_summaryFunction=glbMdlMetricSummaryFn,
#                             model_metric=glbMdlMetricSummary,
#                             model_metric_maximize=glbMdlMetricMaximize)

# Simplify a model
# fit_df <- glbObsFit; glb_mdl <- step(<complex>_mdl)

# Non-caret models
#     rpart_area_mdl <- rpart(reformulate("Area", response=glb_rsp_var), 
#                                data=glbObsFit, #method="class", 
#                                control=rpart.control(cp=0.12),
#                            parms=list(loss=glbMdlMetric_terms))
#     print("rpart_sel_wlm_mdl"); prp(rpart_sel_wlm_mdl)
# 

print(glb_models_df)
```

```
##                                                              id
## MFO###myMFO_classfr                         MFO###myMFO_classfr
## Random###myrandom_classfr             Random###myrandom_classfr
## Max.cor.Y.rcv.1X1###glmnet           Max.cor.Y.rcv.1X1###glmnet
## Max.cor.Y##rcv#rpart                       Max.cor.Y##rcv#rpart
## Max.cor.Y.Time.Poly##rcv#glmnet Max.cor.Y.Time.Poly##rcv#glmnet
## Max.cor.Y.Time.Lag##rcv#glmnet   Max.cor.Y.Time.Lag##rcv#glmnet
## Interact.High.cor.Y##rcv#glmnet Interact.High.cor.Y##rcv#glmnet
## Low.cor.X##rcv#glmnet                     Low.cor.X##rcv#glmnet
## All.X##rcv#glmnet                             All.X##rcv#glmnet
##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              feats
## MFO###myMFO_classfr                                                                                                                                                                                                                                                                                                                                                                                                                                                                         .rnorm
## Random###myrandom_classfr                                                                                                                                                                                                                                                                                                                                                                                                                                                                   .rnorm
## Max.cor.Y.rcv.1X1###glmnet                                                                                                                                                                                                                                                                                                                                                                                                                                            WordCount.root2,NDSS.my.fctr
## Max.cor.Y##rcv#rpart                                                                                                                                                                                                                                                                                                                                                                                                                                                  WordCount.root2,NDSS.my.fctr
## Max.cor.Y.Time.Poly##rcv#glmnet                                                                                                                                                                                                                                                                                                WordCount.root2,NDSS.my.fctr,PubDate.day.minutes.poly.1,PubDate.day.minutes.poly.2,PubDate.day.minutes.poly.3,PubDate.day.minutes.poly.4,PubDate.day.minutes.poly.5
## Max.cor.Y.Time.Lag##rcv#glmnet                                                                                                                                                                                                                                                                                                                                  WordCount.root2,NDSS.my.fctr,PubDate.last2.log1p,PubDate.last4.log1p,PubDate.last8.log1p,PubDate.last16.log1p,PubDate.last32.log1p
## Interact.High.cor.Y##rcv#glmnet                                                                                                                                                                                                                                                                                     WordCount.root2,NDSS.my.fctr,WordCount.root2:WordCount.root2,WordCount.root2:PubDate.day.minutes.poly.1,WordCount.root2:PubDate.last4.log1p,WordCount.root2:PubDate.month.fctr
## Low.cor.X##rcv#glmnet                                                                                                        WordCount.root2,NDSS.my.fctr,PubDate.day.minutes.poly.1,PubDate.wkend,PubDate.day.minutes.poly.4,PubDate.day.minutes.poly.2,PubDate.last4.log1p,PubDate.last16.log1p,PubDate.day.minutes.poly.3,PubDate.month.fctr,.rnorm,PubDate.last32.log1p,PubDate.date.fctr,PubDate.second.fctr,PubDate.minute.fctr,PubDate.wkday.fctr,WordCount.nexp,PubDate.day.minutes.poly.5
## All.X##rcv#glmnet               WordCount.root2,WordCount.log1p,NDSS.my.fctr,PubDate.day.minutes.poly.1,PubDate.hour.fctr,PubDate.wkend,PubDate.day.minutes.poly.4,PubDate.day.minutes.poly.2,PubDate.last4.log1p,PubDate.last2.log1p,PubDate.last8.log1p,PubDate.last16.log1p,PubDate.day.minutes.poly.3,PubDate.month.fctr,PubDate.juliandate,.rnorm,PubDate.last32.log1p,PubDate.date.fctr,PubDate.second.fctr,PubDate.minute.fctr,PubDate.wkday.fctr,WordCount.nexp,PubDate.day.minutes.poly.5
##                                 max.nTuningRuns min.elapsedtime.everything
## MFO###myMFO_classfr                           0                      0.302
## Random###myrandom_classfr                     0                      0.309
## Max.cor.Y.rcv.1X1###glmnet                    0                      1.034
## Max.cor.Y##rcv#rpart                          5                      3.029
## Max.cor.Y.Time.Poly##rcv#glmnet              25                      4.889
## Max.cor.Y.Time.Lag##rcv#glmnet                5                      4.225
## Interact.High.cor.Y##rcv#glmnet              25                      5.157
## Low.cor.X##rcv#glmnet                        25                      6.823
## All.X##rcv#glmnet                            25                      9.521
##                                 min.elapsedtime.final max.AUCpROC.fit
## MFO###myMFO_classfr                             0.004       0.5000000
## Random###myrandom_classfr                       0.002       0.4990604
## Max.cor.Y.rcv.1X1###glmnet                      0.282       0.8790544
## Max.cor.Y##rcv#rpart                            0.077       0.8709432
## Max.cor.Y.Time.Poly##rcv#glmnet                 0.316       0.8748550
## Max.cor.Y.Time.Lag##rcv#glmnet                  0.318       0.8419280
## Interact.High.cor.Y##rcv#glmnet                 0.327       0.8776419
## Low.cor.X##rcv#glmnet                           0.515       0.8751806
## All.X##rcv#glmnet                               0.671       0.8734975
##                                 max.Sens.fit max.Spec.fit max.AUCROCR.fit
## MFO###myMFO_classfr                1.0000000    0.0000000       0.5000000
## Random###myrandom_classfr          0.8312611    0.1668598       0.4972757
## Max.cor.Y.rcv.1X1###glmnet         0.9632073    0.7949015       0.9608594
## Max.cor.Y##rcv#rpart               0.9632073    0.7786790       0.8746354
## Max.cor.Y.Time.Poly##rcv#glmnet    0.9652372    0.7844728       0.9565422
## Max.cor.Y.Time.Lag##rcv#glmnet     0.9781781    0.7056779       0.9581563
## Interact.High.cor.Y##rcv#glmnet    0.9626998    0.7925840       0.9625372
## Low.cor.X##rcv#glmnet              0.9647298    0.7856315       0.9622982
## All.X##rcv#glmnet                  0.9659985    0.7809965       0.9528396
##                                 opt.prob.threshold.fit max.f.score.fit
## MFO###myMFO_classfr                                0.1       0.3045703
## Random###myrandom_classfr                          0.1       0.3045703
## Max.cor.Y.rcv.1X1###glmnet                         0.5       0.8099174
## Max.cor.Y##rcv#rpart                               0.6       0.8000000
## Max.cor.Y.Time.Poly##rcv#glmnet                    0.4       0.8099174
## Max.cor.Y.Time.Lag##rcv#glmnet                     0.3       0.8103347
## Interact.High.cor.Y##rcv#glmnet                    0.4       0.8084359
## Low.cor.X##rcv#glmnet                              0.2       0.8118701
## All.X##rcv#glmnet                                  0.4       0.8101715
##                                 max.Accuracy.fit max.AccuracyLower.fit
## MFO###myMFO_classfr                    0.1796420             0.1688795
## Random###myrandom_classfr              0.1796420             0.1688795
## Max.cor.Y.rcv.1X1###glmnet             0.9329725             0.9255302
## Max.cor.Y##rcv#rpart                   0.9296422             0.9224771
## Max.cor.Y.Time.Poly##rcv#glmnet        0.9322790             0.9255302
## Max.cor.Y.Time.Lag##rcv#glmnet         0.9279084             0.9253120
## Interact.High.cor.Y##rcv#glmnet        0.9315850             0.9244394
## Low.cor.X##rcv#glmnet                  0.9323486             0.9224771
## All.X##rcv#glmnet                      0.9320707             0.9257484
##                                 max.AccuracyUpper.fit max.Kappa.fit
## MFO###myMFO_classfr                         0.1907952     0.0000000
## Random###myrandom_classfr                   0.1907952     0.0000000
## Max.cor.Y.rcv.1X1###glmnet                  0.9398832     0.7692476
## Max.cor.Y##rcv#rpart                        0.9371115     0.7515134
## Max.cor.Y.Time.Poly##rcv#glmnet             0.9398832     0.7643948
## Max.cor.Y.Time.Lag##rcv#glmnet              0.9396854     0.7350030
## Interact.High.cor.Y##rcv#glmnet             0.9388938     0.7641040
## Low.cor.X##rcv#glmnet                       0.9371115     0.7651127
## All.X##rcv#glmnet                           0.9400811     0.7630270
##                                 max.AUCpROC.OOB max.Sens.OOB max.Spec.OOB
## MFO###myMFO_classfr                   0.5000000    1.0000000    0.0000000
## Random###myrandom_classfr             0.5125675    0.8077437    0.2173913
## Max.cor.Y.rcv.1X1###glmnet            0.5962443    0.9098798    0.2826087
## Max.cor.Y##rcv#rpart                  0.5870523    0.9045394    0.2695652
## Max.cor.Y.Time.Poly##rcv#glmnet       0.5965780    0.9105474    0.2826087
## Max.cor.Y.Time.Lag##rcv#glmnet        0.5925379    0.9372497    0.2478261
## Interact.High.cor.Y##rcv#glmnet       0.6009259    0.9105474    0.2913043
## Low.cor.X##rcv#glmnet                 0.6015934    0.9118825    0.2913043
## All.X##rcv#glmnet                     0.5950717    0.9118825    0.2782609
##                                 max.AUCROCR.OOB opt.prob.threshold.OOB
## MFO###myMFO_classfr                   0.5000000                    0.1
## Random###myrandom_classfr             0.4857956                    0.1
## Max.cor.Y.rcv.1X1###glmnet            0.8116126                    0.1
## Max.cor.Y##rcv#rpart                  0.5892132                    0.6
## Max.cor.Y.Time.Poly##rcv#glmnet       0.8049472                    0.1
## Max.cor.Y.Time.Lag##rcv#glmnet        0.8117635                    0.1
## Interact.High.cor.Y##rcv#glmnet       0.8140971                    0.1
## Low.cor.X##rcv#glmnet                 0.8160359                    0.1
## All.X##rcv#glmnet                     0.7995951                    0.1
##                                 max.f.score.OOB max.Accuracy.OOB
## MFO###myMFO_classfr                   0.2349336        0.1331019
## Random###myrandom_classfr             0.2349336        0.1331019
## Max.cor.Y.rcv.1X1###glmnet            0.4405405        0.7604167
## Max.cor.Y##rcv#rpart                  0.2850575        0.8200231
## Max.cor.Y.Time.Poly##rcv#glmnet       0.4441261        0.7754630
## Max.cor.Y.Time.Lag##rcv#glmnet        0.4062196        0.6464120
## Interact.High.cor.Y##rcv#glmnet       0.4398340        0.7656250
## Low.cor.X##rcv#glmnet                 0.4461942        0.7557870
## All.X##rcv#glmnet                     0.4230769        0.7743056
##                                 max.AccuracyLower.OOB
## MFO###myMFO_classfr                         0.1174298
## Random###myrandom_classfr                   0.1174298
## Max.cor.Y.rcv.1X1###glmnet                  0.7395703
## Max.cor.Y##rcv#rpart                        0.8010821
## Max.cor.Y.Time.Poly##rcv#glmnet             0.7550404
## Max.cor.Y.Time.Lag##rcv#glmnet              0.6233476
## Interact.High.cor.Y##rcv#glmnet             0.7449213
## Low.cor.X##rcv#glmnet                       0.7348172
## All.X##rcv#glmnet                           0.7538491
##                                 max.AccuracyUpper.OOB max.Kappa.OOB
## MFO###myMFO_classfr                         0.1500310     0.0000000
## Random###myrandom_classfr                   0.1500310     0.0000000
## Max.cor.Y.rcv.1X1###glmnet                  0.7803749     0.3148374
## Max.cor.Y##rcv#rpart                        0.8378705     0.1825002
## Max.cor.Y.Time.Poly##rcv#glmnet             0.7949457     0.3233542
## Max.cor.Y.Time.Lag##rcv#glmnet              0.6689781     0.2515036
## Interact.High.cor.Y##rcv#glmnet             0.7854227     0.3156027
## Low.cor.X##rcv#glmnet                       0.7758846     0.3197713
## All.X##rcv#glmnet                           0.7938262     0.3001637
##                                 max.AccuracySD.fit max.KappaSD.fit
## MFO###myMFO_classfr                             NA              NA
## Random###myrandom_classfr                       NA              NA
## Max.cor.Y.rcv.1X1###glmnet                      NA              NA
## Max.cor.Y##rcv#rpart                   0.005069520      0.01910910
## Max.cor.Y.Time.Poly##rcv#glmnet        0.005312780      0.01857466
## Max.cor.Y.Time.Lag##rcv#glmnet         0.004559889      0.02005103
## Interact.High.cor.Y##rcv#glmnet        0.005250654      0.01810996
## Low.cor.X##rcv#glmnet                  0.005383249      0.01869090
## All.X##rcv#glmnet                      0.005666511      0.02009517
```

```r
rm(ret_lst)
fit.models_1_chunk_df <- 
    myadd_chunk(fit.models_1_chunk_df, "fit.models_1_end", major.inc = FALSE,
                label.minor = "teardown")
```

```
##                  label step_major step_minor label_minor     bgn     end
## 4 fit.models_1_preProc          1          3     preProc 174.144 174.218
## 5     fit.models_1_end          1          4    teardown 174.218      NA
##   elapsed
## 4   0.074
## 5      NA
```

```r
glb_chunks_df <- myadd_chunk(glb_chunks_df, "fit.models", major.inc = FALSE)
```

```
##         label step_major step_minor label_minor     bgn     end elapsed
## 11 fit.models          6          1           1 150.660 174.228  23.568
## 12 fit.models          6          2           2 174.229      NA      NA
```


```r
fit.models_2_chunk_df <- 
    myadd_chunk(NULL, "fit.models_2_bgn", label.minor = "setup")
```

```
##              label step_major step_minor label_minor     bgn end elapsed
## 1 fit.models_2_bgn          1          0       setup 175.993  NA      NA
```

```r
plt_models_df <- glb_models_df[, -grep("SD|Upper|Lower", names(glb_models_df))]
for (var in grep("^min.", names(plt_models_df), value=TRUE)) {
    plt_models_df[, sub("min.", "inv.", var)] <- 
        #ifelse(all(is.na(tmp <- plt_models_df[, var])), NA, 1.0 / tmp)
        1.0 / plt_models_df[, var]
    plt_models_df <- plt_models_df[ , -grep(var, names(plt_models_df))]
}
print(plt_models_df)
```

```
##                                                              id
## MFO###myMFO_classfr                         MFO###myMFO_classfr
## Random###myrandom_classfr             Random###myrandom_classfr
## Max.cor.Y.rcv.1X1###glmnet           Max.cor.Y.rcv.1X1###glmnet
## Max.cor.Y##rcv#rpart                       Max.cor.Y##rcv#rpart
## Max.cor.Y.Time.Poly##rcv#glmnet Max.cor.Y.Time.Poly##rcv#glmnet
## Max.cor.Y.Time.Lag##rcv#glmnet   Max.cor.Y.Time.Lag##rcv#glmnet
## Interact.High.cor.Y##rcv#glmnet Interact.High.cor.Y##rcv#glmnet
## Low.cor.X##rcv#glmnet                     Low.cor.X##rcv#glmnet
## All.X##rcv#glmnet                             All.X##rcv#glmnet
##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              feats
## MFO###myMFO_classfr                                                                                                                                                                                                                                                                                                                                                                                                                                                                         .rnorm
## Random###myrandom_classfr                                                                                                                                                                                                                                                                                                                                                                                                                                                                   .rnorm
## Max.cor.Y.rcv.1X1###glmnet                                                                                                                                                                                                                                                                                                                                                                                                                                            WordCount.root2,NDSS.my.fctr
## Max.cor.Y##rcv#rpart                                                                                                                                                                                                                                                                                                                                                                                                                                                  WordCount.root2,NDSS.my.fctr
## Max.cor.Y.Time.Poly##rcv#glmnet                                                                                                                                                                                                                                                                                                WordCount.root2,NDSS.my.fctr,PubDate.day.minutes.poly.1,PubDate.day.minutes.poly.2,PubDate.day.minutes.poly.3,PubDate.day.minutes.poly.4,PubDate.day.minutes.poly.5
## Max.cor.Y.Time.Lag##rcv#glmnet                                                                                                                                                                                                                                                                                                                                  WordCount.root2,NDSS.my.fctr,PubDate.last2.log1p,PubDate.last4.log1p,PubDate.last8.log1p,PubDate.last16.log1p,PubDate.last32.log1p
## Interact.High.cor.Y##rcv#glmnet                                                                                                                                                                                                                                                                                     WordCount.root2,NDSS.my.fctr,WordCount.root2:WordCount.root2,WordCount.root2:PubDate.day.minutes.poly.1,WordCount.root2:PubDate.last4.log1p,WordCount.root2:PubDate.month.fctr
## Low.cor.X##rcv#glmnet                                                                                                        WordCount.root2,NDSS.my.fctr,PubDate.day.minutes.poly.1,PubDate.wkend,PubDate.day.minutes.poly.4,PubDate.day.minutes.poly.2,PubDate.last4.log1p,PubDate.last16.log1p,PubDate.day.minutes.poly.3,PubDate.month.fctr,.rnorm,PubDate.last32.log1p,PubDate.date.fctr,PubDate.second.fctr,PubDate.minute.fctr,PubDate.wkday.fctr,WordCount.nexp,PubDate.day.minutes.poly.5
## All.X##rcv#glmnet               WordCount.root2,WordCount.log1p,NDSS.my.fctr,PubDate.day.minutes.poly.1,PubDate.hour.fctr,PubDate.wkend,PubDate.day.minutes.poly.4,PubDate.day.minutes.poly.2,PubDate.last4.log1p,PubDate.last2.log1p,PubDate.last8.log1p,PubDate.last16.log1p,PubDate.day.minutes.poly.3,PubDate.month.fctr,PubDate.juliandate,.rnorm,PubDate.last32.log1p,PubDate.date.fctr,PubDate.second.fctr,PubDate.minute.fctr,PubDate.wkday.fctr,WordCount.nexp,PubDate.day.minutes.poly.5
##                                 max.nTuningRuns max.AUCpROC.fit
## MFO###myMFO_classfr                           0       0.5000000
## Random###myrandom_classfr                     0       0.4990604
## Max.cor.Y.rcv.1X1###glmnet                    0       0.8790544
## Max.cor.Y##rcv#rpart                          5       0.8709432
## Max.cor.Y.Time.Poly##rcv#glmnet              25       0.8748550
## Max.cor.Y.Time.Lag##rcv#glmnet                5       0.8419280
## Interact.High.cor.Y##rcv#glmnet              25       0.8776419
## Low.cor.X##rcv#glmnet                        25       0.8751806
## All.X##rcv#glmnet                            25       0.8734975
##                                 max.Sens.fit max.Spec.fit max.AUCROCR.fit
## MFO###myMFO_classfr                1.0000000    0.0000000       0.5000000
## Random###myrandom_classfr          0.8312611    0.1668598       0.4972757
## Max.cor.Y.rcv.1X1###glmnet         0.9632073    0.7949015       0.9608594
## Max.cor.Y##rcv#rpart               0.9632073    0.7786790       0.8746354
## Max.cor.Y.Time.Poly##rcv#glmnet    0.9652372    0.7844728       0.9565422
## Max.cor.Y.Time.Lag##rcv#glmnet     0.9781781    0.7056779       0.9581563
## Interact.High.cor.Y##rcv#glmnet    0.9626998    0.7925840       0.9625372
## Low.cor.X##rcv#glmnet              0.9647298    0.7856315       0.9622982
## All.X##rcv#glmnet                  0.9659985    0.7809965       0.9528396
##                                 opt.prob.threshold.fit max.f.score.fit
## MFO###myMFO_classfr                                0.1       0.3045703
## Random###myrandom_classfr                          0.1       0.3045703
## Max.cor.Y.rcv.1X1###glmnet                         0.5       0.8099174
## Max.cor.Y##rcv#rpart                               0.6       0.8000000
## Max.cor.Y.Time.Poly##rcv#glmnet                    0.4       0.8099174
## Max.cor.Y.Time.Lag##rcv#glmnet                     0.3       0.8103347
## Interact.High.cor.Y##rcv#glmnet                    0.4       0.8084359
## Low.cor.X##rcv#glmnet                              0.2       0.8118701
## All.X##rcv#glmnet                                  0.4       0.8101715
##                                 max.Accuracy.fit max.Kappa.fit
## MFO###myMFO_classfr                    0.1796420     0.0000000
## Random###myrandom_classfr              0.1796420     0.0000000
## Max.cor.Y.rcv.1X1###glmnet             0.9329725     0.7692476
## Max.cor.Y##rcv#rpart                   0.9296422     0.7515134
## Max.cor.Y.Time.Poly##rcv#glmnet        0.9322790     0.7643948
## Max.cor.Y.Time.Lag##rcv#glmnet         0.9279084     0.7350030
## Interact.High.cor.Y##rcv#glmnet        0.9315850     0.7641040
## Low.cor.X##rcv#glmnet                  0.9323486     0.7651127
## All.X##rcv#glmnet                      0.9320707     0.7630270
##                                 max.AUCpROC.OOB max.Sens.OOB max.Spec.OOB
## MFO###myMFO_classfr                   0.5000000    1.0000000    0.0000000
## Random###myrandom_classfr             0.5125675    0.8077437    0.2173913
## Max.cor.Y.rcv.1X1###glmnet            0.5962443    0.9098798    0.2826087
## Max.cor.Y##rcv#rpart                  0.5870523    0.9045394    0.2695652
## Max.cor.Y.Time.Poly##rcv#glmnet       0.5965780    0.9105474    0.2826087
## Max.cor.Y.Time.Lag##rcv#glmnet        0.5925379    0.9372497    0.2478261
## Interact.High.cor.Y##rcv#glmnet       0.6009259    0.9105474    0.2913043
## Low.cor.X##rcv#glmnet                 0.6015934    0.9118825    0.2913043
## All.X##rcv#glmnet                     0.5950717    0.9118825    0.2782609
##                                 max.AUCROCR.OOB opt.prob.threshold.OOB
## MFO###myMFO_classfr                   0.5000000                    0.1
## Random###myrandom_classfr             0.4857956                    0.1
## Max.cor.Y.rcv.1X1###glmnet            0.8116126                    0.1
## Max.cor.Y##rcv#rpart                  0.5892132                    0.6
## Max.cor.Y.Time.Poly##rcv#glmnet       0.8049472                    0.1
## Max.cor.Y.Time.Lag##rcv#glmnet        0.8117635                    0.1
## Interact.High.cor.Y##rcv#glmnet       0.8140971                    0.1
## Low.cor.X##rcv#glmnet                 0.8160359                    0.1
## All.X##rcv#glmnet                     0.7995951                    0.1
##                                 max.f.score.OOB max.Accuracy.OOB
## MFO###myMFO_classfr                   0.2349336        0.1331019
## Random###myrandom_classfr             0.2349336        0.1331019
## Max.cor.Y.rcv.1X1###glmnet            0.4405405        0.7604167
## Max.cor.Y##rcv#rpart                  0.2850575        0.8200231
## Max.cor.Y.Time.Poly##rcv#glmnet       0.4441261        0.7754630
## Max.cor.Y.Time.Lag##rcv#glmnet        0.4062196        0.6464120
## Interact.High.cor.Y##rcv#glmnet       0.4398340        0.7656250
## Low.cor.X##rcv#glmnet                 0.4461942        0.7557870
## All.X##rcv#glmnet                     0.4230769        0.7743056
##                                 max.Kappa.OOB inv.elapsedtime.everything
## MFO###myMFO_classfr                 0.0000000                  3.3112583
## Random###myrandom_classfr           0.0000000                  3.2362460
## Max.cor.Y.rcv.1X1###glmnet          0.3148374                  0.9671180
## Max.cor.Y##rcv#rpart                0.1825002                  0.3301420
## Max.cor.Y.Time.Poly##rcv#glmnet     0.3233542                  0.2045408
## Max.cor.Y.Time.Lag##rcv#glmnet      0.2515036                  0.2366864
## Interact.High.cor.Y##rcv#glmnet     0.3156027                  0.1939112
## Low.cor.X##rcv#glmnet               0.3197713                  0.1465631
## All.X##rcv#glmnet                   0.3001637                  0.1050310
##                                 inv.elapsedtime.final
## MFO###myMFO_classfr                        250.000000
## Random###myrandom_classfr                  500.000000
## Max.cor.Y.rcv.1X1###glmnet                   3.546099
## Max.cor.Y##rcv#rpart                        12.987013
## Max.cor.Y.Time.Poly##rcv#glmnet              3.164557
## Max.cor.Y.Time.Lag##rcv#glmnet               3.144654
## Interact.High.cor.Y##rcv#glmnet              3.058104
## Low.cor.X##rcv#glmnet                        1.941748
## All.X##rcv#glmnet                            1.490313
```

```r
print(myplot_radar(radar_inp_df=plt_models_df))
```

```
## Warning: The shape palette can deal with a maximum of 6 discrete values
## because more than 6 becomes difficult to discriminate; you have 9.
## Consider specifying shapes manually if you must have them.
```

```
## Warning: Removed 60 rows containing missing values (geom_point).
```

```
## Warning: The shape palette can deal with a maximum of 6 discrete values
## because more than 6 becomes difficult to discriminate; you have 9.
## Consider specifying shapes manually if you must have them.
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_2-1.png) 

```r
# print(myplot_radar(radar_inp_df=subset(plt_models_df, 
#         !(mdl_id %in% grep("random|MFO", plt_models_df$id, value=TRUE)))))

# Compute CI for <metric>SD
glb_models_df <- mutate(glb_models_df, 
                max.df = ifelse(max.nTuningRuns > 1, max.nTuningRuns - 1, NA),
                min.sd2ci.scaler = ifelse(is.na(max.df), NA, qt(0.975, max.df)))
for (var in grep("SD", names(glb_models_df), value=TRUE)) {
    # Does CI alredy exist ?
    var_components <- unlist(strsplit(var, "SD"))
    varActul <- paste0(var_components[1],          var_components[2])
    varUpper <- paste0(var_components[1], "Upper", var_components[2])
    varLower <- paste0(var_components[1], "Lower", var_components[2])
    if (varUpper %in% names(glb_models_df)) {
        warning(varUpper, " already exists in glb_models_df")
        # Assuming Lower also exists
        next
    }    
    print(sprintf("var:%s", var))
    # CI is dependent on sample size in t distribution; df=n-1
    glb_models_df[, varUpper] <- glb_models_df[, varActul] + 
        glb_models_df[, "min.sd2ci.scaler"] * glb_models_df[, var]
    glb_models_df[, varLower] <- glb_models_df[, varActul] - 
        glb_models_df[, "min.sd2ci.scaler"] * glb_models_df[, var]
}
```

```
## Warning: max.AccuracyUpper.fit already exists in glb_models_df
```

```
## [1] "var:max.KappaSD.fit"
```

```r
# Plot metrics with CI
plt_models_df <- glb_models_df[, "id", FALSE]
pltCI_models_df <- glb_models_df[, "id", FALSE]
for (var in grep("Upper", names(glb_models_df), value=TRUE)) {
    var_components <- unlist(strsplit(var, "Upper"))
    col_name <- unlist(paste(var_components, collapse=""))
    plt_models_df[, col_name] <- glb_models_df[, col_name]
    for (name in paste0(var_components[1], c("Upper", "Lower"), var_components[2]))
        pltCI_models_df[, name] <- glb_models_df[, name]
}

build_statsCI_data <- function(plt_models_df) {
    mltd_models_df <- melt(plt_models_df, id.vars="id")
    mltd_models_df$data <- sapply(1:nrow(mltd_models_df), 
        function(row_ix) tail(unlist(strsplit(as.character(
            mltd_models_df[row_ix, "variable"]), "[.]")), 1))
    mltd_models_df$label <- sapply(1:nrow(mltd_models_df), 
        function(row_ix) head(unlist(strsplit(as.character(
            mltd_models_df[row_ix, "variable"]), 
            paste0(".", mltd_models_df[row_ix, "data"]))), 1))
    #print(mltd_models_df)
    
    return(mltd_models_df)
}
mltd_models_df <- build_statsCI_data(plt_models_df)

mltdCI_models_df <- melt(pltCI_models_df, id.vars="id")
for (row_ix in 1:nrow(mltdCI_models_df)) {
    for (type in c("Upper", "Lower")) {
        if (length(var_components <- unlist(strsplit(
                as.character(mltdCI_models_df[row_ix, "variable"]), type))) > 1) {
            #print(sprintf("row_ix:%d; type:%s; ", row_ix, type))
            mltdCI_models_df[row_ix, "label"] <- var_components[1]
            mltdCI_models_df[row_ix, "data"] <- 
                unlist(strsplit(var_components[2], "[.]"))[2]
            mltdCI_models_df[row_ix, "type"] <- type
            break
        }
    }    
}
wideCI_models_df <- reshape(subset(mltdCI_models_df, select=-variable), 
                            timevar="type", 
        idvar=setdiff(names(mltdCI_models_df), c("type", "value", "variable")), 
                            direction="wide")
#print(wideCI_models_df)
mrgdCI_models_df <- merge(wideCI_models_df, mltd_models_df, all.x=TRUE)
#print(mrgdCI_models_df)

# Merge stats back in if CIs don't exist
goback_vars <- c()
for (var in unique(mltd_models_df$label)) {
    for (type in unique(mltd_models_df$data)) {
        var_type <- paste0(var, ".", type)
        # if this data is already present, next
        if (var_type %in% unique(paste(mltd_models_df$label, mltd_models_df$data,
                                       sep=".")))
            next
        #print(sprintf("var_type:%s", var_type))
        goback_vars <- c(goback_vars, var_type)
    }
}

if (length(goback_vars) > 0) {
    mltd_goback_df <- build_statsCI_data(glb_models_df[, c("id", goback_vars)])
    mltd_models_df <- rbind(mltd_models_df, mltd_goback_df)
}

# mltd_models_df <- merge(mltd_models_df, glb_models_df[, c("id", "model_method")], 
#                         all.x=TRUE)

png(paste0(glb_out_pfx, "models_bar.png"), width=480*3, height=480*2)
#print(gp <- myplot_bar(mltd_models_df, "id", "value", colorcol_name="model_method") + 
print(gp <- myplot_bar(df=mltd_models_df, xcol_name="id", ycol_names="value") + 
        geom_errorbar(data=mrgdCI_models_df, 
            mapping=aes(x=mdl_id, ymax=value.Upper, ymin=value.Lower), width=0.5) + 
          facet_grid(label ~ data, scales="free") + 
          theme(axis.text.x = element_text(angle = 90,vjust = 0.5)))
```

```
## Warning: Removed 3 rows containing missing values (geom_errorbar).
```

```r
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
print(gp)
```

```
## Warning: Removed 3 rows containing missing values (geom_errorbar).
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_2-2.png) 

```r
dsp_models_cols <- c("id", 
                    glbMdlMetricsEval[glbMdlMetricsEval %in% names(glb_models_df)],
                    grep("opt.", names(glb_models_df), fixed = TRUE, value = TRUE)) 
# if (glb_is_classification && glb_is_binomial) 
#     dsp_models_cols <- c(dsp_models_cols, "opt.prob.threshold.OOB")
print(dsp_models_df <- orderBy(get_model_sel_frmla(), glb_models_df)[, dsp_models_cols])
```

```
##                                                              id
## Max.cor.Y##rcv#rpart                       Max.cor.Y##rcv#rpart
## Max.cor.Y.Time.Poly##rcv#glmnet Max.cor.Y.Time.Poly##rcv#glmnet
## All.X##rcv#glmnet                             All.X##rcv#glmnet
## Interact.High.cor.Y##rcv#glmnet Interact.High.cor.Y##rcv#glmnet
## Max.cor.Y.rcv.1X1###glmnet           Max.cor.Y.rcv.1X1###glmnet
## Low.cor.X##rcv#glmnet                     Low.cor.X##rcv#glmnet
## Max.cor.Y.Time.Lag##rcv#glmnet   Max.cor.Y.Time.Lag##rcv#glmnet
## MFO###myMFO_classfr                         MFO###myMFO_classfr
## Random###myrandom_classfr             Random###myrandom_classfr
##                                 max.Accuracy.OOB max.AUCROCR.OOB
## Max.cor.Y##rcv#rpart                   0.8200231       0.5892132
## Max.cor.Y.Time.Poly##rcv#glmnet        0.7754630       0.8049472
## All.X##rcv#glmnet                      0.7743056       0.7995951
## Interact.High.cor.Y##rcv#glmnet        0.7656250       0.8140971
## Max.cor.Y.rcv.1X1###glmnet             0.7604167       0.8116126
## Low.cor.X##rcv#glmnet                  0.7557870       0.8160359
## Max.cor.Y.Time.Lag##rcv#glmnet         0.6464120       0.8117635
## MFO###myMFO_classfr                    0.1331019       0.5000000
## Random###myrandom_classfr              0.1331019       0.4857956
##                                 max.AUCpROC.OOB max.Accuracy.fit
## Max.cor.Y##rcv#rpart                  0.5870523        0.9296422
## Max.cor.Y.Time.Poly##rcv#glmnet       0.5965780        0.9322790
## All.X##rcv#glmnet                     0.5950717        0.9320707
## Interact.High.cor.Y##rcv#glmnet       0.6009259        0.9315850
## Max.cor.Y.rcv.1X1###glmnet            0.5962443        0.9329725
## Low.cor.X##rcv#glmnet                 0.6015934        0.9323486
## Max.cor.Y.Time.Lag##rcv#glmnet        0.5925379        0.9279084
## MFO###myMFO_classfr                   0.5000000        0.1796420
## Random###myrandom_classfr             0.5125675        0.1796420
##                                 opt.prob.threshold.fit
## Max.cor.Y##rcv#rpart                               0.6
## Max.cor.Y.Time.Poly##rcv#glmnet                    0.4
## All.X##rcv#glmnet                                  0.4
## Interact.High.cor.Y##rcv#glmnet                    0.4
## Max.cor.Y.rcv.1X1###glmnet                         0.5
## Low.cor.X##rcv#glmnet                              0.2
## Max.cor.Y.Time.Lag##rcv#glmnet                     0.3
## MFO###myMFO_classfr                                0.1
## Random###myrandom_classfr                          0.1
##                                 opt.prob.threshold.OOB
## Max.cor.Y##rcv#rpart                               0.6
## Max.cor.Y.Time.Poly##rcv#glmnet                    0.1
## All.X##rcv#glmnet                                  0.1
## Interact.High.cor.Y##rcv#glmnet                    0.1
## Max.cor.Y.rcv.1X1###glmnet                         0.1
## Low.cor.X##rcv#glmnet                              0.1
## Max.cor.Y.Time.Lag##rcv#glmnet                     0.1
## MFO###myMFO_classfr                                0.1
## Random###myrandom_classfr                          0.1
```

```r
print(myplot_radar(radar_inp_df = dsp_models_df))
```

```
## Warning: The shape palette can deal with a maximum of 6 discrete values
## because more than 6 becomes difficult to discriminate; you have 9.
## Consider specifying shapes manually if you must have them.
```

```
## Warning: Removed 21 rows containing missing values (geom_point).
```

```
## Warning: The shape palette can deal with a maximum of 6 discrete values
## because more than 6 becomes difficult to discriminate; you have 9.
## Consider specifying shapes manually if you must have them.
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_2-3.png) 

```r
print("Metrics used for model selection:"); print(get_model_sel_frmla())
```

```
## [1] "Metrics used for model selection:"
```

```
## ~-max.Accuracy.OOB - max.AUCROCR.OOB - max.AUCpROC.OOB - max.Accuracy.fit - 
##     opt.prob.threshold.OOB
## <environment: 0x7fe3e517d070>
```

```r
print(sprintf("Best model id: %s", dsp_models_df[1, "id"]))
```

```
## [1] "Best model id: Max.cor.Y##rcv#rpart"
```

```r
glb_get_predictions <- function(df, mdl_id, rsp_var, prob_threshold_def=NULL, verbose=FALSE) {
    mdl <- glb_models_lst[[mdl_id]]

    clmnNames <- mygetPredictIds(rsp_var, mdl_id)
    predct_var_name <- clmnNames$value        
    predct_prob_var_name <- clmnNames$prob
    predct_accurate_var_name <- clmnNames$is.acc
    predct_error_var_name <- clmnNames$err
    predct_erabs_var_name <- clmnNames$err.abs

    if (glb_is_regression) {
        df[, predct_var_name] <- predict(mdl, newdata=df, type="raw")
        if (verbose) print(myplot_scatter(df, glb_rsp_var, predct_var_name) + 
                  facet_wrap(reformulate(glbFeatsCategory), scales = "free") + 
                  stat_smooth(method="glm"))

        df[, predct_error_var_name] <- df[, predct_var_name] - df[, glb_rsp_var]
        if (verbose) print(myplot_scatter(df, predct_var_name, predct_error_var_name) + 
                  #facet_wrap(reformulate(glbFeatsCategory), scales = "free") + 
                  stat_smooth(method="auto"))
        if (verbose) print(myplot_scatter(df, glb_rsp_var, predct_error_var_name) + 
                  #facet_wrap(reformulate(glbFeatsCategory), scales = "free") + 
                  stat_smooth(method="glm"))
        
        df[, predct_erabs_var_name] <- abs(df[, predct_error_var_name])
        if (verbose) print(head(orderBy(reformulate(c("-", predct_erabs_var_name)), df)))
        
        df[, predct_accurate_var_name] <- (df[, glb_rsp_var] == df[, predct_var_name])
    }

    if (glb_is_classification && glb_is_binomial) {
        prob_threshold <- glb_models_df[glb_models_df$id == mdl_id, 
                                        "opt.prob.threshold.OOB"]
        if (is.null(prob_threshold) || is.na(prob_threshold)) {
            warning("Using default probability threshold: ", prob_threshold_def)
            if (is.null(prob_threshold <- prob_threshold_def))
                stop("Default probability threshold is NULL")
        }
        
        df[, predct_prob_var_name] <- predict(mdl, newdata = df, type = "prob")[, 2]
        df[, predct_var_name] <- 
        		factor(levels(df[, glb_rsp_var])[
    				(df[, predct_prob_var_name] >=
    					prob_threshold) * 1 + 1], levels(df[, glb_rsp_var]))
    
#         if (verbose) print(myplot_scatter(df, glb_rsp_var, predct_var_name) + 
#                   facet_wrap(reformulate(glbFeatsCategory), scales = "free") + 
#                   stat_smooth(method="glm"))

        df[, predct_error_var_name] <- df[, predct_var_name] != df[, glb_rsp_var]
#         if (verbose) print(myplot_scatter(df, predct_var_name, predct_error_var_name) + 
#                   #facet_wrap(reformulate(glbFeatsCategory), scales = "free") + 
#                   stat_smooth(method="auto"))
#         if (verbose) print(myplot_scatter(df, glb_rsp_var, predct_error_var_name) + 
#                   #facet_wrap(reformulate(glbFeatsCategory), scales = "free") + 
#                   stat_smooth(method="glm"))
        
        # if prediction is a TP (true +ve), measure distance from 1.0
        tp <- which((df[, predct_var_name] == df[, glb_rsp_var]) &
                    (df[, predct_var_name] == levels(df[, glb_rsp_var])[2]))
        df[tp, predct_erabs_var_name] <- abs(1 - df[tp, predct_prob_var_name])
        #rowIx <- which.max(df[tp, predct_erabs_var_name]); df[tp, c(glb_id_var, glb_rsp_var, predct_var_name, predct_prob_var_name, predct_erabs_var_name)][rowIx, ]
        
        # if prediction is a TN (true -ve), measure distance from 0.0
        tn <- which((df[, predct_var_name] == df[, glb_rsp_var]) &
                    (df[, predct_var_name] == levels(df[, glb_rsp_var])[1]))
        df[tn, predct_erabs_var_name] <- abs(0 - df[tn, predct_prob_var_name])
        #rowIx <- which.max(df[tn, predct_erabs_var_name]); df[tn, c(glb_id_var, glb_rsp_var, predct_var_name, predct_prob_var_name, predct_erabs_var_name)][rowIx, ]
        
        # if prediction is a FP (flse +ve), measure distance from 0.0
        fp <- which((df[, predct_var_name] != df[, glb_rsp_var]) &
                    (df[, predct_var_name] == levels(df[, glb_rsp_var])[2]))
        df[fp, predct_erabs_var_name] <- abs(0 - df[fp, predct_prob_var_name])
        #rowIx <- which.max(df[fp, predct_erabs_var_name]); df[fp, c(glb_id_var, glb_rsp_var, predct_var_name, predct_prob_var_name, predct_erabs_var_name)][rowIx, ]
        
        # if prediction is a FN (flse -ve), measure distance from 1.0
        fn <- which((df[, predct_var_name] != df[, glb_rsp_var]) &
                    (df[, predct_var_name] == levels(df[, glb_rsp_var])[1]))
        df[fn, predct_erabs_var_name] <- abs(1 - df[fn, predct_prob_var_name])
        #rowIx <- which.max(df[fn, predct_erabs_var_name]); df[fn, c(glb_id_var, glb_rsp_var, predct_var_name, predct_prob_var_name, predct_erabs_var_name)][rowIx, ]

        
        if (verbose) print(head(orderBy(reformulate(c("-", predct_erabs_var_name)), df)))
        
        df[, predct_accurate_var_name] <- (df[, glb_rsp_var] == df[, predct_var_name])
    }    
    
    if (glb_is_classification && !glb_is_binomial) {
        df[, predct_var_name] <- predict(mdl, newdata = df, type = "raw")
        df[, paste0(predct_var_name, ".prob")] <- 
            predict(mdl, newdata = df, type = "prob")
        stop("Multinomial prediction error calculation needs to be implemented...")
    }

    return(df)
}    

#stop(here"); glb2Sav(); glbObsAll <- savObsAll; glbObsTrn <- savObsTrn; glbObsFit <- savObsFit; glbObsOOB <- savObsOOB; sav_models_df <- glb_models_df; glb_models_df <- sav_models_df; glb_featsimp_df <- sav_featsimp_df    

myget_category_stats <- function(obs_df, mdl_id, label) {
    require(dplyr)
    require(lazyeval)
    
    predct_var_name <- mygetPredictIds(glb_rsp_var, mdl_id)$value        
    predct_error_var_name <- mygetPredictIds(glb_rsp_var, mdl_id)$err.abs
    
    if (!predct_var_name %in% names(obs_df))
        obs_df <- glb_get_predictions(obs_df, mdl_id, glb_rsp_var)
    
    tmp_obs_df <- obs_df[, c(glbFeatsCategory, glb_rsp_var, 
                             predct_var_name, predct_error_var_name)]
#     tmp_obs_df <- obs_df %>%
#         dplyr::select_(glbFeatsCategory, glb_rsp_var, predct_var_name, predct_error_var_name) 
    #dplyr::rename(startprice.log10.predict.RFE.X.glmnet.err=error_abs_OOB)
    names(tmp_obs_df)[length(names(tmp_obs_df))] <- paste0("err.abs.", label)
    
    ret_ctgry_df <- tmp_obs_df %>%
        dplyr::group_by_(glbFeatsCategory) %>%
        dplyr::summarise_(#interp(~sum(abs(var)), var=as.name(glb_rsp_var)), 
            interp(~sum(var), var=as.name(paste0("err.abs.", label))), 
            interp(~mean(var), var=as.name(paste0("err.abs.", label))),
            interp(~n()))
    names(ret_ctgry_df) <- c(glbFeatsCategory, 
                             #paste0(glb_rsp_var, ".abs.", label, ".sum"),
                             paste0("err.abs.", label, ".sum"),         
                             paste0("err.abs.", label, ".mean"), 
                             paste0(".n.", label))
    ret_ctgry_df <- dplyr::ungroup(ret_ctgry_df)
    #colSums(ret_ctgry_df[, -grep(glbFeatsCategory, names(ret_ctgry_df))])
    
    return(ret_ctgry_df)    
}
#print(colSums((ctgry_df <- myget_category_stats(obs_df=glbObsFit, mdl_id="", label="fit"))[, -grep(glbFeatsCategory, names(ctgry_df))]))

if (!is.null(glb_mdl_ensemble)) {
    fit.models_2_chunk_df <- myadd_chunk(fit.models_2_chunk_df, 
                            paste0("fit.models_2_", mdl_id_pfx), major.inc = TRUE, 
                                                label.minor = "ensemble")
    
    mdl_id_pfx <- "Ensemble"

    if (#(glb_is_regression) | 
        ((glb_is_classification) & (!glb_is_binomial)))
        stop("Ensemble models not implemented yet for multinomial classification")
    
    mygetEnsembleAutoMdlIds <- function() {
        tmp_models_df <- orderBy(get_model_sel_frmla(), glb_models_df)
        row.names(tmp_models_df) <- tmp_models_df$id
        mdl_threshold_pos <- 
            min(which(grepl("MFO|Random|Baseline", tmp_models_df$id))) - 1
        mdlIds <- tmp_models_df$id[1:mdl_threshold_pos]
        return(mdlIds[!grepl("Ensemble", mdlIds)])
    }
    
    if (glb_mdl_ensemble == "auto") {
        glb_mdl_ensemble <- mygetEnsembleAutoMdlIds()
        mdl_id_pfx <- paste0(mdl_id_pfx, ".auto")        
    } else if (grepl("^%<d-%", glb_mdl_ensemble)) {
        glb_mdl_ensemble <- eval(parse(text =
                        str_trim(unlist(strsplit(glb_mdl_ensemble, "%<d-%"))[2])))
    }
    
    for (mdl_id in glb_mdl_ensemble) {
        if (!(mdl_id %in% names(glb_models_lst))) {
            warning("Model ", mdl_id, " in glb_model_ensemble not found !")
            next
        }
        glbObsFit <- glb_get_predictions(df = glbObsFit, mdl_id, glb_rsp_var)
        glbObsOOB <- glb_get_predictions(df = glbObsOOB, mdl_id, glb_rsp_var)
    }
    
#mdl_id_pfx <- "Ensemble.RFE"; mdlId <- paste0(mdl_id_pfx, ".glmnet")
#glb_mdl_ensemble <- gsub(mygetPredictIds$value, "", grep("RFE\\.X\\.(?!Interact)", row.names(glb_featsimp_df), perl = TRUE, value = TRUE), fixed = TRUE)
#varImp(glb_models_lst[[mdlId]])
    
#cor_df <- data.frame(cor=cor(glbObsFit[, glb_rsp_var], glbObsFit[, paste(mygetPredictIds$value, glb_mdl_ensemble)], use="pairwise.complete.obs"))
#glbObsFit <- glb_get_predictions(df=glbObsFit, "Ensemble.glmnet", glb_rsp_var);print(colSums((ctgry_df <- myget_category_stats(obs_df=glbObsFit, mdl_id="Ensemble.glmnet", label="fit"))[, -grep(glbFeatsCategory, names(ctgry_df))]))
    
    ### bid0_sp
    #  Better than MFO; models.n=28; min.RMSE.fit=0.0521233; err.abs.fit.sum=7.3631895
    #  old: Top x from auto; models.n= 5; min.RMSE.fit=0.06311047; err.abs.fit.sum=9.5937080
    #  RFE only ;       models.n=16; min.RMSE.fit=0.05148588; err.abs.fit.sum=7.2875091
    #  RFE subset only ;models.n= 5; min.RMSE.fit=0.06040702; err.abs.fit.sum=9.059088
    #  RFE subset only ;models.n= 9; min.RMSE.fit=0.05933167; err.abs.fit.sum=8.7421288
    #  RFE subset only ;models.n=15; min.RMSE.fit=0.0584607; err.abs.fit.sum=8.5902066
    #  RFE subset only ;models.n=17; min.RMSE.fit=0.05496899; err.abs.fit.sum=8.0170431
    #  RFE subset only ;models.n=18; min.RMSE.fit=0.05441577; err.abs.fit.sum=7.837223
    #  RFE subset only ;models.n=16; min.RMSE.fit=0.05441577; err.abs.fit.sum=7.837223
    ### bid0_sp
    ### bid1_sp
    # "auto"; err.abs.fit.sum=76.699774; min.RMSE.fit=0.2186429
    # "RFE.X.*"; err.abs.fit.sum=; min.RMSE.fit=0.221114
    ### bid1_sp

    indep_vars <- paste(mygetPredictIds(glb_rsp_var)$value, glb_mdl_ensemble, sep = "")
    if (glb_is_classification)
        indep_vars <- paste(indep_vars, ".prob", sep = "")
    # Some models in glb_mdl_ensemble might not be fitted e.g. RFE.X.Interact
    indep_vars <- intersect(indep_vars, names(glbObsFit))
    
#     indep_vars <- grep(mygetPredictIds(glb_rsp_var)$value, names(glbObsFit), fixed=TRUE, value=TRUE)
#     if (glb_is_regression)
#         indep_vars <- indep_vars[!grepl("(err\\.abs|accurate)$", indep_vars)]
#     if (glb_is_classification && glb_is_binomial)
#         indep_vars <- grep("prob$", indep_vars, value=TRUE) else
#         indep_vars <- indep_vars[!grepl("err$", indep_vars)]

    #rfe_fit_ens_results <- myrun_rfe(glbObsFit, indep_vars)
    
    for (method in c("glm", "glmnet")) {
        for (trainControlMethod in 
             c("boot", "boot632", "cv", "repeatedcv"
               #, "LOOCV" # tuneLength * nrow(fitDF)
               , "LGOCV", "adaptive_cv"
               #, "adaptive_boot"  #error: adaptive$min should be less than 3 
               #, "adaptive_LGOCV" #error: adaptive$min should be less than 3 
               )) {
            #sav_models_df <- glb_models_df; all.equal(sav_models_df, glb_models_df)
            #glb_models_df <- sav_models_df; print(glb_models_df$id)
                
            if ((method == "glm") && (trainControlMethod != "repeatedcv"))
                # glm used only to identify outliers
                next
            
            ret_lst <- myfit_mdl(
                mdl_specs_lst = myinit_mdl_specs_lst(mdl_specs_lst = list(
                    id.prefix = paste0(mdl_id_pfx, ".", trainControlMethod), 
                    type = glb_model_type, tune.df = NULL,
                    trainControl.method = trainControlMethod,
                    trainControl.number = glb_rcv_n_folds,
                    trainControl.repeats = glb_rcv_n_repeats,
                    trainControl.classProbs = glb_is_classification,
                    trainControl.summaryFunction = glbMdlMetricSummaryFn,
                    train.metric = glbMdlMetricSummary, 
                    train.maximize = glbMdlMetricMaximize,    
                    train.method = method)),
                indep_vars = indep_vars, rsp_var = glb_rsp_var, 
                fit_df = glbObsFit, OOB_df = glbObsOOB)
        }
    }
    dsp_models_df <- get_dsp_models_df()
}

if (is.null(glb_sel_mdl_id)) 
    glb_sel_mdl_id <- dsp_models_df[1, "id"] else 
    print(sprintf("User specified selection: %s", glb_sel_mdl_id))   
```

```
## [1] "User specified selection: All.X##rcv#glmnet"
```

```r
myprint_mdl(glb_sel_mdl <- glb_models_lst[[glb_sel_mdl_id]])
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_2-4.png) 

```
##             Length Class      Mode     
## a0            89   -none-     numeric  
## beta        5073   dgCMatrix  S4       
## df            89   -none-     numeric  
## dim            2   -none-     numeric  
## lambda        89   -none-     numeric  
## dev.ratio     89   -none-     numeric  
## nulldev        1   -none-     numeric  
## npasses        1   -none-     numeric  
## jerr           1   -none-     numeric  
## offset         1   -none-     logical  
## classnames     2   -none-     character
## call           5   -none-     call     
## nobs           1   -none-     numeric  
## lambdaOpt      1   -none-     numeric  
## xNames        57   -none-     character
## problemType    1   -none-     character
## tuneValue      2   data.frame list     
## obsLevels      2   -none-     character
## [1] "min lambda > lambdaOpt:"
##                     (Intercept)       NDSS.my.fctr#Opnn#RmFrDbt 
##                     -4.15778549                     -0.29068581 
##    NDSS.my.fctr#Opnn#ThPblcEdtr NDSS.my.fctrBsnss#Crsswrds/Gms# 
##                      2.70584441                      3.68584256 
##          NDSS.my.fctrOpEd#Opnn#          NDSS.my.fctrScnc#Hlth# 
##                      4.00374690                      3.09723359 
##         NDSS.my.fctrStyls#U.S.#             NDSS.my.fctrTStyl## 
##                      2.85084048                     -0.13236677 
##      PubDate.day.minutes.poly.1                 WordCount.log1p 
##                      5.59078306                      0.08550402 
##                 WordCount.root2 
##                      0.04496835 
## [1] "max lambda < lambdaOpt:"
##                     (Intercept)       NDSS.my.fctr#Opnn#RmFrDbt 
##                    -4.284190646                    -0.447707141 
##    NDSS.my.fctr#Opnn#ThPblcEdtr NDSS.my.fctrBsnss#Crsswrds/Gms# 
##                     2.827765103                     3.751516261 
##      NDSS.my.fctrBsnss#Tchnlgy#          NDSS.my.fctrOpEd#Opnn# 
##                     0.081535690                     4.070511321 
##          NDSS.my.fctrScnc#Hlth#         NDSS.my.fctrStyls#U.S.# 
##                     3.163208431                     2.914061199 
##             NDSS.my.fctrTStyl##      PubDate.day.minutes.poly.1 
##                    -0.187607543                     7.038308760 
##                   PubDate.wkend                 WordCount.log1p 
##                     0.008146344                     0.096036181 
##                 WordCount.root2 
##                     0.046479071
```

```
## [1] TRUE
```

```r
# From here to save(), this should all be in one function
#   these are executed in the same seq twice more:
#       fit.data.training & predict.data.new chunks
print(sprintf("%s fit prediction diagnostics:", glb_sel_mdl_id))
```

```
## [1] "All.X##rcv#glmnet fit prediction diagnostics:"
```

```r
glbObsFit <- glb_get_predictions(df = glbObsFit, mdl_id = glb_sel_mdl_id, 
                                 rsp_var = glb_rsp_var)
print(sprintf("%s OOB prediction diagnostics:", glb_sel_mdl_id))
```

```
## [1] "All.X##rcv#glmnet OOB prediction diagnostics:"
```

```r
glbObsOOB <- glb_get_predictions(df = glbObsOOB, mdl_id = glb_sel_mdl_id, 
                                     rsp_var = glb_rsp_var)

glb_featsimp_df <- 
    myget_feats_importance(mdl=glb_sel_mdl, featsimp_df=NULL)
glb_featsimp_df[, paste0(glb_sel_mdl_id, ".imp")] <- glb_featsimp_df$imp
#mdl_id <-"RFE.X.glmnet"; glb_featsimp_df <- myget_feats_importance(glb_models_lst[[mdl_id]], glb_featsimp_df); glb_featsimp_df[, paste0(mdl_id, ".imp")] <- glb_featsimp_df$imp; print(glb_featsimp_df)
#print(head(sbst_featsimp_df <- subset(glb_featsimp_df, is.na(RFE.X.glmnet.imp) | (abs(RFE.X.YeoJohnson.glmnet.imp - RFE.X.glmnet.imp) > 0.0001), select=-imp)))
#print(orderBy(~ -cor.y.abs, subset(glb_feats_df, id %in% c(row.names(sbst_featsimp_df), "startprice.dcm1.is9", "D.weight.post.stop.sum"))))
print(glb_featsimp_df)
```

```
##                                            imp All.X##rcv#glmnet.imp
## PubDate.day.minutes.poly.1          100.000000            100.000000
## NDSS.my.fctrOpEd#Opnn#               68.562523             68.562523
## NDSS.my.fctrBsnss#Crsswrds/Gms#      63.559761             63.559761
## NDSS.my.fctrScnc#Hlth#               54.307908             54.307908
## NDSS.my.fctrStyls#U.S.#              50.421526             50.421526
## NDSS.my.fctr#Opnn#ThPblcEdtr         48.417943             48.417943
## WordCount.log1p                       6.701272              6.701272
## WordCount.root2                       6.021636              6.021636
## NDSS.my.fctrBsnss#Tchnlgy#            5.690755              5.690755
## PubDate.wkend                         5.345889              5.345889
## .rnorm                                5.307609              5.307609
## NDSS.my.fctr#Mltmd#                   5.307609              5.307609
## NDSS.my.fctr#U.S.#Edctn               5.307609              5.307609
## NDSS.my.fctrBsnss#BsnssDy#Dlbk        5.307609              5.307609
## NDSS.my.fctrBsnss#BsnssDy#SmllBsnss   5.307609              5.307609
## NDSS.my.fctrCltr##                    5.307609              5.307609
## NDSS.my.fctrCltr#Arts#                5.307609              5.307609
## NDSS.my.fctrFrgn#Wrld#                5.307609              5.307609
## NDSS.my.fctrFrgn#Wrld#AsPcfc          5.307609              5.307609
## NDSS.my.fctrMtr#N.Y./Rgn#             5.307609              5.307609
## NDSS.my.fctrStyls##Fshn               5.307609              5.307609
## NDSS.my.fctrTrvl#Trvl#                5.307609              5.307609
## NDSS.my.fctrmyOthr                    5.307609              5.307609
## PubDate.date.fctr(7,13]               5.307609              5.307609
## PubDate.date.fctr(13,19]              5.307609              5.307609
## PubDate.date.fctr(19,25]              5.307609              5.307609
## PubDate.date.fctr(25,31]              5.307609              5.307609
## PubDate.day.minutes.poly.2            5.307609              5.307609
## PubDate.day.minutes.poly.3            5.307609              5.307609
## PubDate.day.minutes.poly.4            5.307609              5.307609
## PubDate.day.minutes.poly.5            5.307609              5.307609
## PubDate.hour.fctr(7.67,15.3]          5.307609              5.307609
## PubDate.hour.fctr(15.3,23]            5.307609              5.307609
## PubDate.juliandate                    5.307609              5.307609
## PubDate.last16.log1p                  5.307609              5.307609
## PubDate.last2.log1p                   5.307609              5.307609
## PubDate.last32.log1p                  5.307609              5.307609
## PubDate.last4.log1p                   5.307609              5.307609
## PubDate.last8.log1p                   5.307609              5.307609
## PubDate.minute.fctr(14.8,29.5]        5.307609              5.307609
## PubDate.minute.fctr(29.5,44.2]        5.307609              5.307609
## PubDate.minute.fctr(44.2,59.1]        5.307609              5.307609
## PubDate.month.fctr10                  5.307609              5.307609
## PubDate.month.fctr11                  5.307609              5.307609
## PubDate.month.fctr12                  5.307609              5.307609
## PubDate.second.fctr(14.8,29.5]        5.307609              5.307609
## PubDate.second.fctr(29.5,44.2]        5.307609              5.307609
## PubDate.second.fctr(44.2,59.1]        5.307609              5.307609
## PubDate.wkday.fctr1                   5.307609              5.307609
## PubDate.wkday.fctr2                   5.307609              5.307609
## PubDate.wkday.fctr3                   5.307609              5.307609
## PubDate.wkday.fctr4                   5.307609              5.307609
## PubDate.wkday.fctr5                   5.307609              5.307609
## PubDate.wkday.fctr6                   5.307609              5.307609
## WordCount.nexp                        5.307609              5.307609
## NDSS.my.fctrTStyl##                   2.967145              2.967145
## NDSS.my.fctr#Opnn#RmFrDbt             0.000000              0.000000
```

```r
# Used again in fit.data.training & predict.data.new chunks
glb_analytics_diag_plots <- function(obs_df, mdl_id, prob_threshold=NULL) {
    if (!is.null(featsimp_df <- glb_featsimp_df)) {
        featsimp_df$feat <- gsub("`(.*?)`", "\\1", row.names(featsimp_df))    
        featsimp_df$feat.interact <- gsub("(.*?):(.*)", "\\2", featsimp_df$feat)
        featsimp_df$feat <- gsub("(.*?):(.*)", "\\1", featsimp_df$feat)    
        featsimp_df$feat.interact <- 
            ifelse(featsimp_df$feat.interact == featsimp_df$feat, 
                                            NA, featsimp_df$feat.interact)
        featsimp_df$feat <- 
            gsub("(.*?)\\.fctr(.*)", "\\1\\.fctr", featsimp_df$feat)
        featsimp_df$feat.interact <- 
            gsub("(.*?)\\.fctr(.*)", "\\1\\.fctr", featsimp_df$feat.interact) 
        featsimp_df <- orderBy(~ -imp.max, 
            summaryBy(imp ~ feat + feat.interact, data=featsimp_df,
                      FUN=max))    
        #rex_str=":(.*)"; txt_vctr=tail(featsimp_df$feat); ret_lst <- regexec(rex_str, txt_vctr); ret_lst <- regmatches(txt_vctr, ret_lst); ret_vctr <- sapply(1:length(ret_lst), function(pos_ix) ifelse(length(ret_lst[[pos_ix]]) > 0, ret_lst[[pos_ix]], "")); print(ret_vctr <- ret_vctr[ret_vctr != ""])    
        
        featsimp_df <- subset(featsimp_df, !is.na(imp.max))
        if (nrow(featsimp_df) > 5) {
            warning("Limiting important feature scatter plots to 5 out of ",
                    nrow(featsimp_df))
            featsimp_df <- head(featsimp_df, 5)
        }
        
    #     if (!all(is.na(featsimp_df$feat.interact)))
    #         stop("not implemented yet")
        rsp_var_out <- mygetPredictIds(glb_rsp_var, mdl_id)$value
        for (var in featsimp_df$feat) {
            plot_df <- melt(obs_df, id.vars = var, 
                            measure.vars = c(glb_rsp_var, rsp_var_out))
    
            print(myplot_scatter(plot_df, var, "value", colorcol_name = "variable",
                                facet_colcol_name = "variable", jitter = TRUE) + 
                          guides(color = FALSE))
        }
    }
    
    if (glb_is_regression) {
        if (is.null(featsimp_df) || (nrow(featsimp_df) == 0))
            warning("No important features in glb_fin_mdl") else
            print(myplot_prediction_regression(df=obs_df, 
                        feat_x=ifelse(nrow(featsimp_df) > 1, featsimp_df$feat[2],
                                      ".rownames"), 
                                               feat_y=featsimp_df$feat[1],
                        rsp_var=glb_rsp_var, rsp_var_out=rsp_var_out,
                        id_vars=glb_id_var)
    #               + facet_wrap(reformulate(featsimp_df$feat[2])) # if [1 or 2] is a factor
    #               + geom_point(aes_string(color="<col_name>.fctr")) #  to color the plot
                  )
    }    
    
    if (glb_is_classification) {
        if (is.null(featsimp_df) || (nrow(featsimp_df) == 0))
            warning("No features in selected model are statistically important")
        else print(myplot_prediction_classification(df = obs_df, 
                                feat_x = ifelse(nrow(featsimp_df) > 1, 
                                                featsimp_df$feat[2], ".rownames"),
                                               feat_y = featsimp_df$feat[1],
                                                rsp_var = glb_rsp_var, 
                                                rsp_var_out = rsp_var_out, 
                                                id_vars = glb_id_var,
                                                prob_threshold = prob_threshold))
    }    
}

if (glb_is_classification && glb_is_binomial)
    glb_analytics_diag_plots(obs_df = glbObsOOB, mdl_id = glb_sel_mdl_id, 
            prob_threshold = glb_models_df[glb_models_df$id == glb_sel_mdl_id, 
                                           "opt.prob.threshold.OOB"]) else
    glb_analytics_diag_plots(obs_df = glbObsOOB, mdl_id = glb_sel_mdl_id)                  
```

```
## Warning in glb_analytics_diag_plots(obs_df = glbObsOOB, mdl_id =
## glb_sel_mdl_id, : Limiting important feature scatter plots to 5 out of 23
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_2-5.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_2-6.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_2-7.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_2-8.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_2-9.png) 

```
## [1] "Min/Max Boundaries: "
##   UniqueID Pplr.fctr Pplr.fctr.All.X..rcv.glmnet.prob
## 1     2555         N                       0.01260243
## 2      302         N                       0.10198151
##   Pplr.fctr.All.X..rcv.glmnet Pplr.fctr.All.X..rcv.glmnet.err
## 1                           N                           FALSE
## 2                           Y                            TRUE
##   Pplr.fctr.All.X..rcv.glmnet.err.abs Pplr.fctr.All.X..rcv.glmnet.is.acc
## 1                          0.01260243                               TRUE
## 2                          0.10198151                              FALSE
##   Pplr.fctr.All.X..rcv.glmnet.accurate Pplr.fctr.All.X..rcv.glmnet.error
## 1                                 TRUE                       0.000000000
## 2                                FALSE                       0.001981509
##   .label
## 1   2555
## 2    302
## [1] "Inaccurate: "
##   UniqueID Pplr.fctr Pplr.fctr.All.X..rcv.glmnet.prob
## 1     4020         Y                       0.02959503
## 2     4775         Y                       0.03504702
## 3     6354         Y                       0.03540579
## 4     4745         Y                       0.03706549
## 5      172         Y                       0.03819629
## 6     2283         Y                       0.03988159
##   Pplr.fctr.All.X..rcv.glmnet Pplr.fctr.All.X..rcv.glmnet.err
## 1                           N                            TRUE
## 2                           N                            TRUE
## 3                           N                            TRUE
## 4                           N                            TRUE
## 5                           N                            TRUE
## 6                           N                            TRUE
##   Pplr.fctr.All.X..rcv.glmnet.err.abs Pplr.fctr.All.X..rcv.glmnet.is.acc
## 1                           0.9704050                              FALSE
## 2                           0.9649530                              FALSE
## 3                           0.9645942                              FALSE
## 4                           0.9629345                              FALSE
## 5                           0.9618037                              FALSE
## 6                           0.9601184                              FALSE
##   Pplr.fctr.All.X..rcv.glmnet.accurate Pplr.fctr.All.X..rcv.glmnet.error
## 1                                FALSE                       -0.07040497
## 2                                FALSE                       -0.06495298
## 3                                FALSE                       -0.06459421
## 4                                FALSE                       -0.06293451
## 5                                FALSE                       -0.06180371
## 6                                FALSE                       -0.06011841
##     UniqueID Pplr.fctr Pplr.fctr.All.X..rcv.glmnet.prob
## 94      3334         N                        0.1012627
## 200      982         N                        0.1389648
## 218     1500         N                        0.1565391
## 238     5832         N                        0.2091344
## 242     5441         N                        0.2151956
## 383     1448         N                        0.9039811
##     Pplr.fctr.All.X..rcv.glmnet Pplr.fctr.All.X..rcv.glmnet.err
## 94                            Y                            TRUE
## 200                           Y                            TRUE
## 218                           Y                            TRUE
## 238                           Y                            TRUE
## 242                           Y                            TRUE
## 383                           Y                            TRUE
##     Pplr.fctr.All.X..rcv.glmnet.err.abs Pplr.fctr.All.X..rcv.glmnet.is.acc
## 94                            0.1012627                              FALSE
## 200                           0.1389648                              FALSE
## 218                           0.1565391                              FALSE
## 238                           0.2091344                              FALSE
## 242                           0.2151956                              FALSE
## 383                           0.9039811                              FALSE
##     Pplr.fctr.All.X..rcv.glmnet.accurate Pplr.fctr.All.X..rcv.glmnet.error
## 94                                 FALSE                       0.001262711
## 200                                FALSE                       0.038964799
## 218                                FALSE                       0.056539106
## 238                                FALSE                       0.109134369
## 242                                FALSE                       0.115195645
## 383                                FALSE                       0.803981082
##     UniqueID Pplr.fctr Pplr.fctr.All.X..rcv.glmnet.prob
## 385     4943         N                        0.9166355
## 386      770         N                        0.9216170
## 387      221         N                        0.9266596
## 388     3590         N                        0.9284441
## 389      472         N                        0.9287613
## 390     2995         N                        0.9308453
##     Pplr.fctr.All.X..rcv.glmnet Pplr.fctr.All.X..rcv.glmnet.err
## 385                           Y                            TRUE
## 386                           Y                            TRUE
## 387                           Y                            TRUE
## 388                           Y                            TRUE
## 389                           Y                            TRUE
## 390                           Y                            TRUE
##     Pplr.fctr.All.X..rcv.glmnet.err.abs Pplr.fctr.All.X..rcv.glmnet.is.acc
## 385                           0.9166355                              FALSE
## 386                           0.9216170                              FALSE
## 387                           0.9266596                              FALSE
## 388                           0.9284441                              FALSE
## 389                           0.9287613                              FALSE
## 390                           0.9308453                              FALSE
##     Pplr.fctr.All.X..rcv.glmnet.accurate Pplr.fctr.All.X..rcv.glmnet.error
## 385                                FALSE                         0.8166355
## 386                                FALSE                         0.8216170
## 387                                FALSE                         0.8266596
## 388                                FALSE                         0.8284441
## 389                                FALSE                         0.8287613
## 390                                FALSE                         0.8308453
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_2-10.png) 

```r
if (!is.null(glbFeatsCategory)) {
    glbLvlCategory <- merge(glbLvlCategory, 
            myget_category_stats(obs_df = glbObsFit, mdl_id = glb_sel_mdl_id, 
                                 label = "fit"), 
                            by = glbFeatsCategory, all = TRUE)
    row.names(glbLvlCategory) <- glbLvlCategory[, glbFeatsCategory]
    glbLvlCategory <- merge(glbLvlCategory, 
            myget_category_stats(obs_df = glbObsOOB, mdl_id = glb_sel_mdl_id,
                                 label="OOB"),
                          #by=glbFeatsCategory, all=TRUE) glb_ctgry-df already contains .n.OOB ?
                          all = TRUE)
    row.names(glbLvlCategory) <- glbLvlCategory[, glbFeatsCategory]
    if (any(grepl("OOB", glbMdlMetricsEval)))
        print(orderBy(~-err.abs.OOB.mean, glbLvlCategory)) else
            print(orderBy(~-err.abs.fit.mean, glbLvlCategory))
    print(colSums(glbLvlCategory[, -grep(glbFeatsCategory, names(glbLvlCategory))]))
}
```

```
##                                    NDSS.my.fctr .n.OOB .n.Fit .n.Tst
## OpEd#Opnn#                           OpEd#Opnn#     89    437    164
## Styls#U.S.#                         Styls#U.S.#     50    127     61
## Scnc#Hlth#                           Scnc#Hlth#     48    148     57
## #Opnn#ThPblcEdtr               #Opnn#ThPblcEdtr      4     16     10
## Bsnss#Crsswrds/Gms#         Bsnss#Crsswrds/Gms#     18    105     42
## Bsnss#Tchnlgy#                   Bsnss#Tchnlgy#    126    213    114
## #Opnn#RmFrDbt                     #Opnn#RmFrDbt     20     42     20
## Bsnss#BsnssDy#Dlbk           Bsnss#BsnssDy#Dlbk    323    629    304
## ##                                           ##    371    913    342
## Cltr#Arts#                           Cltr#Arts#    185    490    174
## Mtr#N.Y./Rgn#                     Mtr#N.Y./Rgn#     70    128     67
## Styls##Fshn                         Styls##Fshn     15    104     15
## Bsnss#BsnssDy#SmllBsnss Bsnss#BsnssDy#SmllBsnss     40    100     41
## Frgn#Wrld#AsPcfc               Frgn#Wrld#AsPcfc     53    150     56
## TStyl##                                 TStyl##    101    623    105
## Trvl#Trvl#                           Trvl#Trvl#     34     83     35
## #Mltmd#                                 #Mltmd#     49     92     52
## myOthr                                   myOthr      5     33      5
## #U.S.#Edctn                         #U.S.#Edctn     82    243     89
## Cltr##                                   Cltr##      1     NA     70
## Frgn#Wrld#                           Frgn#Wrld#     44    128     47
##                         .freqRatio.Fit .freqRatio.OOB .freqRatio.Tst
## OpEd#Opnn#                 0.090965862   0.0515046296    0.087700535
## Styls#U.S.#                0.026436303   0.0289351852    0.032620321
## Scnc#Hlth#                 0.030807660   0.0277777778    0.030481283
## #Opnn#ThPblcEdtr           0.003330558   0.0023148148    0.005347594
## Bsnss#Crsswrds/Gms#        0.021856786   0.0104166667    0.022459893
## Bsnss#Tchnlgy#             0.044338052   0.0729166667    0.060962567
## #Opnn#RmFrDbt              0.008742714   0.0115740741    0.010695187
## Bsnss#BsnssDy#Dlbk         0.130932556   0.1869212963    0.162566845
## ##                         0.190049958   0.2146990741    0.182887701
## Cltr#Arts#                 0.101998335   0.1070601852    0.093048128
## Mtr#N.Y./Rgn#              0.026644463   0.0405092593    0.035828877
## Styls##Fshn                0.021648626   0.0086805556    0.008021390
## Bsnss#BsnssDy#SmllBsnss    0.020815987   0.0231481481    0.021925134
## Frgn#Wrld#AsPcfc           0.031223980   0.0306712963    0.029946524
## TStyl##                    0.129683597   0.0584490741    0.056149733
## Trvl#Trvl#                 0.017277269   0.0196759259    0.018716578
## #Mltmd#                    0.019150708   0.0283564815    0.027807487
## myOthr                     0.006869276   0.0028935185    0.002673797
## #U.S.#Edctn                0.050582848   0.0474537037    0.047593583
## Cltr##                              NA   0.0005787037    0.037433155
## Frgn#Wrld#                 0.026644463   0.0254629630    0.025133690
##                         err.abs.fit.sum err.abs.fit.mean .n.fit
## OpEd#Opnn#                   119.836470       0.27422533    437
## Styls#U.S.#                   56.628946       0.44589721    127
## Scnc#Hlth#                    57.015658       0.38524093    148
## #Opnn#ThPblcEdtr               6.584894       0.41155588     16
## Bsnss#Crsswrds/Gms#           27.017444       0.25730899    105
## Bsnss#Tchnlgy#                39.873202       0.18719813    213
## #Opnn#RmFrDbt                  7.034814       0.16749556     42
## Bsnss#BsnssDy#Dlbk            81.819044       0.13007797    629
## ##                           101.017935       0.11064396    913
## Cltr#Arts#                    48.216560       0.09840114    490
## Mtr#N.Y./Rgn#                 15.400336       0.12031512    128
## Styls##Fshn                    6.787573       0.06526513    104
## Bsnss#BsnssDy#SmllBsnss       10.499902       0.10499902    100
## Frgn#Wrld#AsPcfc              13.904763       0.09269842    150
## TStyl##                       33.400369       0.05361215    623
## Trvl#Trvl#                     4.002682       0.04822508     83
## #Mltmd#                        5.892332       0.06404709     92
## myOthr                         2.320443       0.07031647     33
## #U.S.#Edctn                   12.184453       0.05014178    243
## Cltr##                               NA               NA     NA
## Frgn#Wrld#                     5.444890       0.04253820    128
##                         err.abs.OOB.sum err.abs.OOB.mean
## OpEd#Opnn#                  51.18065179       0.57506350
## Styls#U.S.#                 25.99760848       0.51995217
## Scnc#Hlth#                  24.69049636       0.51438534
## #Opnn#ThPblcEdtr             1.99109686       0.49777421
## Bsnss#Crsswrds/Gms#          8.73848753       0.48547153
## Bsnss#Tchnlgy#              25.76388022       0.20447524
## #Opnn#RmFrDbt                3.91776433       0.19588822
## Bsnss#BsnssDy#Dlbk          59.90568728       0.18546652
## ##                          67.73033042       0.18256154
## Cltr#Arts#                  31.26604266       0.16900564
## Mtr#N.Y./Rgn#               11.39362621       0.16276609
## Styls##Fshn                  1.91928194       0.12795213
## Bsnss#BsnssDy#SmllBsnss      4.62878656       0.11571966
## Frgn#Wrld#AsPcfc             4.79443402       0.09046102
## TStyl##                      7.86776348       0.07789865
## Trvl#Trvl#                   2.55254529       0.07507486
## #Mltmd#                      3.51781240       0.07179209
## myOthr                       0.33933658       0.06786732
## #U.S.#Edctn                  5.11442705       0.06237106
## Cltr##                       0.04526133       0.04526133
## Frgn#Wrld#                   1.88150036       0.04276137
##           .n.OOB           .n.Fit           .n.Tst   .freqRatio.Fit 
##      1728.000000               NA      1870.000000               NA 
##   .freqRatio.OOB   .freqRatio.Tst  err.abs.fit.sum err.abs.fit.mean 
##         1.000000         1.000000               NA               NA 
##           .n.fit  err.abs.OOB.sum err.abs.OOB.mean 
##               NA       345.236821         4.469969
```

```r
write.csv(glbObsOOB[, c(glb_id_var, 
                grep(glb_rsp_var, names(glbObsOOB), fixed=TRUE, value=TRUE))], 
    paste0(gsub(".", "_", paste0(glb_out_pfx, glb_sel_mdl_id), fixed=TRUE), 
           "_OOBobs.csv"), row.names=FALSE)

fit.models_2_chunk_df <- 
    myadd_chunk(NULL, "fit.models_2_bgn", label.minor = "teardown")
```

```
##              label step_major step_minor label_minor     bgn end elapsed
## 1 fit.models_2_bgn          1          0    teardown 187.068  NA      NA
```

```r
glb_chunks_df <- myadd_chunk(glb_chunks_df, "fit.models", major.inc=FALSE)
```

```
##         label step_major step_minor label_minor     bgn     end elapsed
## 12 fit.models          6          2           2 174.229 187.078  12.849
## 13 fit.models          6          3           3 187.078      NA      NA
```


```r
# if (sum(is.na(glbObsAll$D.P.http)) > 0)
#         stop("fit.models_3: Why is this happening ?")

#stop(here"); glb2Sav()
sync_glb_obs_df <- function() {
    # Merge or cbind ?
    for (col in setdiff(names(glbObsFit), names(glbObsTrn)))
        glbObsTrn[glbObsTrn$.lcn == "Fit", col] <<- glbObsFit[, col]
    for (col in setdiff(names(glbObsFit), names(glbObsAll)))
        glbObsAll[glbObsAll$.lcn == "Fit", col] <<- glbObsFit[, col]
    if (all(is.na(glbObsNew[, glb_rsp_var])))
        for (col in setdiff(names(glbObsOOB), names(glbObsTrn)))
            glbObsTrn[glbObsTrn$.lcn == "OOB", col] <<- glbObsOOB[, col]
    for (col in setdiff(names(glbObsOOB), names(glbObsAll)))
        glbObsAll[glbObsAll$.lcn == "OOB", col] <<- glbObsOOB[, col]
}
sync_glb_obs_df()
    
print(setdiff(names(glbObsNew), names(glbObsAll)))
```

```
## character(0)
```

```r
if (glb_save_envir)
    save(glb_feats_df, 
         glbObsAll, #glbObsTrn, glbObsFit, glbObsOOB, glbObsNew,
         glb_models_df, dsp_models_df, glb_models_lst, glb_sel_mdl, glb_sel_mdl_id,
         glb_model_type,
        file=paste0(glb_out_pfx, "selmdl_dsk.RData"))
#load(paste0(glb_out_pfx, "selmdl_dsk.RData"))

rm(ret_lst)
```

```
## Warning in rm(ret_lst): object 'ret_lst' not found
```

```r
replay.petrisim(pn=glb_analytics_pn, 
    replay.trans=(glb_analytics_avl_objs <- c(glb_analytics_avl_objs, 
        "model.selected")), flip_coord=TRUE)
```

```
## time	trans	 "bgn " "fit.data.training.all " "predict.data.new " "end " 
## 0.0000 	multiple enabled transitions:  data.training.all data.new model.selected 	firing:  data.training.all 
## 1.0000 	 1 	 2 1 0 0 
## 1.0000 	multiple enabled transitions:  data.training.all data.new model.selected model.final data.training.all.prediction 	firing:  data.new 
## 2.0000 	 2 	 1 1 1 0 
## 2.0000 	multiple enabled transitions:  data.training.all data.new model.selected model.final data.training.all.prediction data.new.prediction 	firing:  model.selected 
## 3.0000 	 3 	 0 2 1 0
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_3-1.png) 

```r
glb_chunks_df <- myadd_chunk(glb_chunks_df, "fit.data.training", major.inc=TRUE)
```

```
##                label step_major step_minor label_minor     bgn     end
## 13        fit.models          6          3           3 187.078 192.947
## 14 fit.data.training          7          0           0 192.948      NA
##    elapsed
## 13   5.869
## 14      NA
```

## Step `7.0: fit data training`

```r
#load(paste0(glb_inp_pfx, "dsk.RData"))

#stop(here"); glb2Sav()
if (!is.null(glb_fin_mdl_id) && (glb_fin_mdl_id %in% names(glb_models_lst))) {
    warning("Final model same as user selected model")
    glb_fin_mdl <- glb_models_lst[[glb_fin_mdl_id]]
} else 
# if (nrow(glbObsFit) + length(glbObsFitOutliers) == nrow(glbObsTrn))
if (!all(is.na(glbObsNew[, glb_rsp_var])))
{    
    warning("Final model same as glb_sel_mdl_id")
    glb_fin_mdl_id <- paste0("Final.", glb_sel_mdl_id)
    glb_fin_mdl <- glb_sel_mdl
    glb_models_lst[[glb_fin_mdl_id]] <- glb_fin_mdl
} else {    
    
#     if (grepl("RFE", glb_sel_mdl_id) || 
#         (!is.null(glb_mdl_ensemble) && grepl("RFE", glb_mdl_ensemble))) {
        indep_vars <- myadjust_interaction_feats(subset(glb_feats_df, 
                                            !nzv & (exclude.as.feat != 1))[, "id"])
        rfe_trn_results <- 
            myrun_rfe(glbObsTrn, indep_vars, glbRFESizes[["Final"]])
        if (!isTRUE(all.equal(sort(predictors(rfe_trn_results)),
                              sort(predictors(rfe_fit_results))))) {
            print("Diffs predictors(rfe_trn_results) vs. predictors(rfe_fit_results):")
            print(setdiff(predictors(rfe_trn_results), predictors(rfe_fit_results)))
            print("Diffs predictors(rfe_fit_results) vs. predictors(rfe_trn_results):")
            print(setdiff(predictors(rfe_fit_results), predictors(rfe_trn_results)))
        }
    # }    

    if (grepl("Ensemble", glb_sel_mdl_id)) {
        # Find which models are relevant
        mdlimp_df <- subset(myget_feats_importance(glb_sel_mdl), imp > 5)
        # Fit selected models on glbObsTrn
        for (mdl_id in gsub(".prob", "", 
gsub(mygetPredictIds(glb_rsp_var)$value, "", row.names(mdlimp_df), fixed = TRUE),
                            fixed = TRUE)) {
            mdl_id_components <- unlist(strsplit(mdl_id, "[.]"))
            mdlIdPfx <- paste0(c(head(mdl_id_components, -1), "Train"), 
                               collapse = ".")
            if (grepl("RFE\\.X\\.", mdlIdPfx)) 
                mdlIndepVars <- myadjust_interaction_feats(myextract_actual_feats(
                    predictors(rfe_trn_results))) else
                mdlIndepVars <- trim(unlist(
            strsplit(glb_models_df[glb_models_df$id == mdl_id, "feats"], "[,]")))
            ret_lst <- 
                myfit_mdl(mdl_specs_lst = myinit_mdl_specs_lst(mdl_specs_lst = list(
                        id.prefix = mdlIdPfx, 
                        type = glb_model_type, tune.df = glbMdlTuneParams,
                        trainControl.method = "repeatedcv",
                        trainControl.number = glb_rcv_n_folds,
                        trainControl.repeats = glb_rcv_n_repeats,
                        trainControl.classProbs = glb_is_classification,
                        trainControl.summaryFunction = glbMdlMetricSummaryFn,
                        train.metric = glbMdlMetricSummary, 
                        train.maximize = glbMdlMetricMaximize,    
                        train.method = tail(mdl_id_components, 1))),
                    indep_vars = mdlIndepVars,
                    rsp_var = glb_rsp_var, 
                    fit_df = glbObsTrn, OOB_df = NULL)
            
            glbObsTrn <- glb_get_predictions(df = glbObsTrn,
                                                mdl_id = tail(glb_models_df$id, 1), 
                                                rsp_var = glb_rsp_var,
                                                prob_threshold_def = 
                    subset(glb_models_df, id == mdl_id)$opt.prob.threshold.OOB)
            glbObsNew <- glb_get_predictions(df = glbObsNew,
                                                mdl_id = tail(glb_models_df$id, 1), 
                                                rsp_var = glb_rsp_var,
                                                prob_threshold_def = 
                    subset(glb_models_df, id == mdl_id)$opt.prob.threshold.OOB)
        }    
    }
    
    # "Final" model
    if ((model_method <- glb_sel_mdl$method) == "custom")
        # get actual method from the mdl_id
        model_method <- tail(unlist(strsplit(glb_sel_mdl_id, "[.]")), 1)
        
    if (grepl("Ensemble", glb_sel_mdl_id)) {
        # Find which models are relevant
        mdlimp_df <- subset(myget_feats_importance(glb_sel_mdl), imp > 5)
        if (glb_is_classification && glb_is_binomial)
            indep_vars_vctr <- gsub("(.*)\\.(.*)\\.prob", "\\1\\.Train\\.\\2\\.prob",
                                    row.names(mdlimp_df)) else
            indep_vars_vctr <- gsub("(.*)\\.(.*)", "\\1\\.Train\\.\\2",
                                    row.names(mdlimp_df))
    } else 
    if (grepl("RFE.X", glb_sel_mdl_id, fixed = TRUE)) {
        indep_vars_vctr <- myextract_actual_feats(predictors(rfe_trn_results))
    } else indep_vars_vctr <- 
                trim(unlist(strsplit(glb_models_df[glb_models_df$id ==
                                                   glb_sel_mdl_id
                                                   , "feats"], "[,]")))
        
    if (!is.null(glb_preproc_methods) &&
        ((match_pos <- regexpr(gsub(".", "\\.", 
                                    paste(glb_preproc_methods, collapse = "|"),
                                   fixed = TRUE), glb_sel_mdl_id)) != -1))
        ths_preProcess <- str_sub(glb_sel_mdl_id, match_pos, 
                                match_pos + attr(match_pos, "match.length") - 1) else
        ths_preProcess <- NULL                                      

    mdl_id_pfx <- ifelse(grepl("Ensemble", glb_sel_mdl_id),
                                   "Final.Ensemble", "Final")
    trnobs_df <- if (is.null(glbObsTrnOutliers[[mdl_id_pfx]])) glbObsTrn else 
        glbObsTrn[!(glbObsTrn[, glb_id_var] %in%
                            glbObsTrnOutliers[[mdl_id_pfx]]), ]
        
    # Force fitting of Final.glm to identify outliers
    #method_vctr <- unique(c("glm", myparseMdlId(glb_sel_mdl_id)$alg))
    # or skip glm for speed
    method_vctr <- myparseMdlId(glb_sel_mdl_id)$alg
    for (method in method_vctr) {
        #source("caret_nominalTrainWorkflow.R")
        
        # glmnet requires at least 2 indep vars
        if ((length(indep_vars_vctr) == 1) && (method %in% "glmnet"))
            next

        ret_lst <- 
            myfit_mdl(mdl_specs_lst = myinit_mdl_specs_lst(mdl_specs_lst = list(
                id.prefix = mdl_id_pfx, 
                type = glb_model_type, trainControl.method = "repeatedcv",
                trainControl.number = glb_rcv_n_folds, 
                trainControl.repeats = glb_rcv_n_repeats,
                trainControl.classProbs = glb_is_classification,
                trainControl.summaryFunction = glbMdlMetricSummaryFn,
trainControl.allowParallel = if (method %in% c("glm", "glmnet")) FALSE else TRUE,
                train.metric = glbMdlMetricSummary, 
                train.maximize = glbMdlMetricMaximize,    
                train.method = method,
                train.preProcess = ths_preProcess)),
                indep_vars = indep_vars_vctr, rsp_var = glb_rsp_var, 
                fit_df = trnobs_df, OOB_df = NULL)
    }
        
    if ((length(method_vctr) == 1) || (method != "glm")) {
        glb_fin_mdl <- glb_models_lst[[length(glb_models_lst)]] 
        glb_fin_mdl_id <- glb_models_df[length(glb_models_lst), "id"]
    }
}
```

```
## +(rfe) fit Fold1.Rep1 size: 55
```

```
## Warning in lda.default(x, grouping, ...): variables are collinear
```

```
## -(rfe) fit Fold1.Rep1 size: 55 
## +(rfe) imp Fold1.Rep1 
## -(rfe) imp Fold1.Rep1 
## +(rfe) fit Fold1.Rep1 size: 32 
## -(rfe) fit Fold1.Rep1 size: 32 
## +(rfe) fit Fold1.Rep1 size: 16 
## -(rfe) fit Fold1.Rep1 size: 16 
## +(rfe) fit Fold1.Rep1 size:  8 
## -(rfe) fit Fold1.Rep1 size:  8 
## +(rfe) fit Fold1.Rep1 size:  4 
## -(rfe) fit Fold1.Rep1 size:  4 
## +(rfe) fit Fold1.Rep1 size:  2 
## -(rfe) fit Fold1.Rep1 size:  2 
## +(rfe) fit Fold2.Rep1 size: 55
```

```
## Warning in lda.default(x, grouping, ...): variables are collinear
```

```
## -(rfe) fit Fold2.Rep1 size: 55 
## +(rfe) imp Fold2.Rep1 
## -(rfe) imp Fold2.Rep1 
## +(rfe) fit Fold2.Rep1 size: 32 
## -(rfe) fit Fold2.Rep1 size: 32 
## +(rfe) fit Fold2.Rep1 size: 16 
## -(rfe) fit Fold2.Rep1 size: 16 
## +(rfe) fit Fold2.Rep1 size:  8 
## -(rfe) fit Fold2.Rep1 size:  8 
## +(rfe) fit Fold2.Rep1 size:  4 
## -(rfe) fit Fold2.Rep1 size:  4 
## +(rfe) fit Fold2.Rep1 size:  2 
## -(rfe) fit Fold2.Rep1 size:  2 
## +(rfe) fit Fold3.Rep1 size: 55
```

```
## Warning in lda.default(x, grouping, ...): variables are collinear
```

```
## -(rfe) fit Fold3.Rep1 size: 55 
## +(rfe) imp Fold3.Rep1 
## -(rfe) imp Fold3.Rep1 
## +(rfe) fit Fold3.Rep1 size: 32 
## -(rfe) fit Fold3.Rep1 size: 32 
## +(rfe) fit Fold3.Rep1 size: 16 
## -(rfe) fit Fold3.Rep1 size: 16 
## +(rfe) fit Fold3.Rep1 size:  8 
## -(rfe) fit Fold3.Rep1 size:  8 
## +(rfe) fit Fold3.Rep1 size:  4 
## -(rfe) fit Fold3.Rep1 size:  4 
## +(rfe) fit Fold3.Rep1 size:  2 
## -(rfe) fit Fold3.Rep1 size:  2 
## +(rfe) fit Fold1.Rep2 size: 55
```

```
## Warning in lda.default(x, grouping, ...): variables are collinear
```

```
## -(rfe) fit Fold1.Rep2 size: 55 
## +(rfe) imp Fold1.Rep2 
## -(rfe) imp Fold1.Rep2 
## +(rfe) fit Fold1.Rep2 size: 32 
## -(rfe) fit Fold1.Rep2 size: 32 
## +(rfe) fit Fold1.Rep2 size: 16 
## -(rfe) fit Fold1.Rep2 size: 16 
## +(rfe) fit Fold1.Rep2 size:  8 
## -(rfe) fit Fold1.Rep2 size:  8 
## +(rfe) fit Fold1.Rep2 size:  4 
## -(rfe) fit Fold1.Rep2 size:  4 
## +(rfe) fit Fold1.Rep2 size:  2 
## -(rfe) fit Fold1.Rep2 size:  2 
## +(rfe) fit Fold2.Rep2 size: 55
```

```
## Warning in lda.default(x, grouping, ...): variables are collinear
```

```
## -(rfe) fit Fold2.Rep2 size: 55 
## +(rfe) imp Fold2.Rep2 
## -(rfe) imp Fold2.Rep2 
## +(rfe) fit Fold2.Rep2 size: 32 
## -(rfe) fit Fold2.Rep2 size: 32 
## +(rfe) fit Fold2.Rep2 size: 16 
## -(rfe) fit Fold2.Rep2 size: 16 
## +(rfe) fit Fold2.Rep2 size:  8 
## -(rfe) fit Fold2.Rep2 size:  8 
## +(rfe) fit Fold2.Rep2 size:  4 
## -(rfe) fit Fold2.Rep2 size:  4 
## +(rfe) fit Fold2.Rep2 size:  2 
## -(rfe) fit Fold2.Rep2 size:  2 
## +(rfe) fit Fold3.Rep2 size: 55
```

```
## Warning in lda.default(x, grouping, ...): variables are collinear
```

```
## -(rfe) fit Fold3.Rep2 size: 55 
## +(rfe) imp Fold3.Rep2 
## -(rfe) imp Fold3.Rep2 
## +(rfe) fit Fold3.Rep2 size: 32 
## -(rfe) fit Fold3.Rep2 size: 32 
## +(rfe) fit Fold3.Rep2 size: 16 
## -(rfe) fit Fold3.Rep2 size: 16 
## +(rfe) fit Fold3.Rep2 size:  8 
## -(rfe) fit Fold3.Rep2 size:  8 
## +(rfe) fit Fold3.Rep2 size:  4 
## -(rfe) fit Fold3.Rep2 size:  4 
## +(rfe) fit Fold3.Rep2 size:  2 
## -(rfe) fit Fold3.Rep2 size:  2 
## +(rfe) fit Fold1.Rep3 size: 55
```

```
## Warning in lda.default(x, grouping, ...): variables are collinear
```

```
## -(rfe) fit Fold1.Rep3 size: 55 
## +(rfe) imp Fold1.Rep3 
## -(rfe) imp Fold1.Rep3 
## +(rfe) fit Fold1.Rep3 size: 32 
## -(rfe) fit Fold1.Rep3 size: 32 
## +(rfe) fit Fold1.Rep3 size: 16 
## -(rfe) fit Fold1.Rep3 size: 16 
## +(rfe) fit Fold1.Rep3 size:  8 
## -(rfe) fit Fold1.Rep3 size:  8 
## +(rfe) fit Fold1.Rep3 size:  4 
## -(rfe) fit Fold1.Rep3 size:  4 
## +(rfe) fit Fold1.Rep3 size:  2 
## -(rfe) fit Fold1.Rep3 size:  2 
## +(rfe) fit Fold2.Rep3 size: 55
```

```
## Warning in lda.default(x, grouping, ...): variables are collinear
```

```
## -(rfe) fit Fold2.Rep3 size: 55 
## +(rfe) imp Fold2.Rep3 
## -(rfe) imp Fold2.Rep3 
## +(rfe) fit Fold2.Rep3 size: 32 
## -(rfe) fit Fold2.Rep3 size: 32 
## +(rfe) fit Fold2.Rep3 size: 16 
## -(rfe) fit Fold2.Rep3 size: 16 
## +(rfe) fit Fold2.Rep3 size:  8 
## -(rfe) fit Fold2.Rep3 size:  8 
## +(rfe) fit Fold2.Rep3 size:  4 
## -(rfe) fit Fold2.Rep3 size:  4 
## +(rfe) fit Fold2.Rep3 size:  2 
## -(rfe) fit Fold2.Rep3 size:  2 
## +(rfe) fit Fold3.Rep3 size: 55
```

```
## Warning in lda.default(x, grouping, ...): variables are collinear
```

```
## -(rfe) fit Fold3.Rep3 size: 55 
## +(rfe) imp Fold3.Rep3 
## -(rfe) imp Fold3.Rep3 
## +(rfe) fit Fold3.Rep3 size: 32 
## -(rfe) fit Fold3.Rep3 size: 32 
## +(rfe) fit Fold3.Rep3 size: 16 
## -(rfe) fit Fold3.Rep3 size: 16 
## +(rfe) fit Fold3.Rep3 size:  8 
## -(rfe) fit Fold3.Rep3 size:  8 
## +(rfe) fit Fold3.Rep3 size:  4 
## -(rfe) fit Fold3.Rep3 size:  4 
## +(rfe) fit Fold3.Rep3 size:  2 
## -(rfe) fit Fold3.Rep3 size:  2
```

```
## Warning in lda.default(x, grouping, ...): variables are collinear
```

```
## 
## Recursive feature selection
## 
## Outer resampling method: Cross-Validated (3 fold, repeated 3 times) 
## 
## Resampling performance over subset size:
## 
##  Variables Accuracy   Kappa AccuracySD KappaSD Selected
##          2   0.8204 0.03718   0.003445 0.01035         
##          4   0.8760 0.44502   0.003122 0.01919         
##          8   0.8738 0.44420   0.003283 0.02061         
##         16   0.9016 0.63787   0.006552 0.02352         
##         32   0.9011 0.63667   0.006474 0.02303         
##         55   0.9030 0.64644   0.006845 0.02410        *
## 
## The top 5 variables (out of 55):
##    WordCount.log1p, WordCount.root2, WordCount.nexp, NDSS.my.fctrOpEd#Opnn#, PubDate.day.minutes.poly.1
## 
##  [1] "WordCount.log1p"                    
##  [2] "WordCount.root2"                    
##  [3] "WordCount.nexp"                     
##  [4] "NDSS.my.fctrOpEd#Opnn#"             
##  [5] "PubDate.day.minutes.poly.1"         
##  [6] "PubDate.day.minutes.poly.4"         
##  [7] "PubDate.hour.fctr(15.3,23]"         
##  [8] "PubDate.last4.log1p"                
##  [9] "PubDate.last2.log1p"                
## [10] "NDSS.my.fctrScnc#Hlth#"             
## [11] "NDSS.my.fctrBsnss#Crsswrds/Gms#"    
## [12] "PubDate.day.minutes.poly.5"         
## [13] "PubDate.last8.log1p"                
## [14] "NDSS.my.fctrStyls#U.S.#"            
## [15] "PubDate.wkend"                      
## [16] "PubDate.last16.log1p"               
## [17] "PubDate.day.minutes.poly.2"         
## [18] "PubDate.juliandate"                 
## [19] "PubDate.wkday.fctr6"                
## [20] "PubDate.month.fctr11"               
## [21] "PubDate.second.fctr(14.8,29.5]"     
## [22] "PubDate.date.fctr(7,13]"            
## [23] ".rnorm"                             
## [24] "PubDate.wkday.fctr1"                
## [25] "PubDate.day.minutes.poly.3"         
## [26] "PubDate.date.fctr(25,31]"           
## [27] "PubDate.hour.fctr(7.67,15.3]"       
## [28] "PubDate.last32.log1p"               
## [29] "PubDate.minute.fctr(14.8,29.5]"     
## [30] "PubDate.month.fctr10"               
## [31] "NDSS.my.fctrBsnss#Tchnlgy#"         
## [32] "PubDate.second.fctr(29.5,44.2]"     
## [33] "NDSS.my.fctrmyOthr"                 
## [34] "PubDate.date.fctr(13,19]"           
## [35] "PubDate.wkday.fctr3"                
## [36] "PubDate.minute.fctr(44.2,59.1]"     
## [37] "PubDate.wkday.fctr4"                
## [38] "PubDate.second.fctr(44.2,59.1]"     
## [39] "NDSS.my.fctr#Opnn#RmFrDbt"          
## [40] "PubDate.date.fctr(19,25]"           
## [41] "NDSS.my.fctrMtr#N.Y./Rgn#"          
## [42] "NDSS.my.fctrBsnss#BsnssDy#SmllBsnss"
## [43] "NDSS.my.fctrTrvl#Trvl#"             
## [44] "NDSS.my.fctrStyls##Fshn"            
## [45] "NDSS.my.fctr#Mltmd#"                
## [46] "PubDate.wkday.fctr2"                
## [47] "NDSS.my.fctrFrgn#Wrld#"             
## [48] "NDSS.my.fctrFrgn#Wrld#AsPcfc"       
## [49] "PubDate.wkday.fctr5"                
## [50] "PubDate.minute.fctr(29.5,44.2]"     
## [51] "NDSS.my.fctr#U.S.#Edctn"            
## [52] "NDSS.my.fctrCltr#Arts#"             
## [53] "NDSS.my.fctrBsnss#BsnssDy#Dlbk"     
## [54] "NDSS.my.fctr##"                     
## [55] "NDSS.my.fctrTStyl##"
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.data.training_0-1.png) 

```
## [1] "fitting model: Final##rcv#glmnet"
## [1] "    indep_vars: WordCount.root2,WordCount.log1p,NDSS.my.fctr,PubDate.day.minutes.poly.1,PubDate.hour.fctr,PubDate.wkend,PubDate.day.minutes.poly.4,PubDate.day.minutes.poly.2,PubDate.last4.log1p,PubDate.last2.log1p,PubDate.last8.log1p,PubDate.last16.log1p,PubDate.day.minutes.poly.3,PubDate.month.fctr,PubDate.juliandate,.rnorm,PubDate.last32.log1p,PubDate.date.fctr,PubDate.second.fctr,PubDate.minute.fctr,PubDate.wkday.fctr,WordCount.nexp,PubDate.day.minutes.poly.5"
## + Fold1.Rep1: alpha=0.100, lambda=0.07781 
## - Fold1.Rep1: alpha=0.100, lambda=0.07781 
## + Fold1.Rep1: alpha=0.325, lambda=0.07781 
## - Fold1.Rep1: alpha=0.325, lambda=0.07781 
## + Fold1.Rep1: alpha=0.550, lambda=0.07781 
## - Fold1.Rep1: alpha=0.550, lambda=0.07781 
## + Fold1.Rep1: alpha=0.775, lambda=0.07781 
## - Fold1.Rep1: alpha=0.775, lambda=0.07781 
## + Fold1.Rep1: alpha=1.000, lambda=0.07781 
## - Fold1.Rep1: alpha=1.000, lambda=0.07781 
## + Fold2.Rep1: alpha=0.100, lambda=0.07781 
## - Fold2.Rep1: alpha=0.100, lambda=0.07781 
## + Fold2.Rep1: alpha=0.325, lambda=0.07781 
## - Fold2.Rep1: alpha=0.325, lambda=0.07781 
## + Fold2.Rep1: alpha=0.550, lambda=0.07781 
## - Fold2.Rep1: alpha=0.550, lambda=0.07781 
## + Fold2.Rep1: alpha=0.775, lambda=0.07781 
## - Fold2.Rep1: alpha=0.775, lambda=0.07781 
## + Fold2.Rep1: alpha=1.000, lambda=0.07781 
## - Fold2.Rep1: alpha=1.000, lambda=0.07781 
## + Fold3.Rep1: alpha=0.100, lambda=0.07781 
## - Fold3.Rep1: alpha=0.100, lambda=0.07781 
## + Fold3.Rep1: alpha=0.325, lambda=0.07781 
## - Fold3.Rep1: alpha=0.325, lambda=0.07781 
## + Fold3.Rep1: alpha=0.550, lambda=0.07781 
## - Fold3.Rep1: alpha=0.550, lambda=0.07781 
## + Fold3.Rep1: alpha=0.775, lambda=0.07781 
## - Fold3.Rep1: alpha=0.775, lambda=0.07781 
## + Fold3.Rep1: alpha=1.000, lambda=0.07781 
## - Fold3.Rep1: alpha=1.000, lambda=0.07781 
## + Fold1.Rep2: alpha=0.100, lambda=0.07781 
## - Fold1.Rep2: alpha=0.100, lambda=0.07781 
## + Fold1.Rep2: alpha=0.325, lambda=0.07781 
## - Fold1.Rep2: alpha=0.325, lambda=0.07781 
## + Fold1.Rep2: alpha=0.550, lambda=0.07781 
## - Fold1.Rep2: alpha=0.550, lambda=0.07781 
## + Fold1.Rep2: alpha=0.775, lambda=0.07781 
## - Fold1.Rep2: alpha=0.775, lambda=0.07781 
## + Fold1.Rep2: alpha=1.000, lambda=0.07781 
## - Fold1.Rep2: alpha=1.000, lambda=0.07781 
## + Fold2.Rep2: alpha=0.100, lambda=0.07781 
## - Fold2.Rep2: alpha=0.100, lambda=0.07781 
## + Fold2.Rep2: alpha=0.325, lambda=0.07781 
## - Fold2.Rep2: alpha=0.325, lambda=0.07781 
## + Fold2.Rep2: alpha=0.550, lambda=0.07781 
## - Fold2.Rep2: alpha=0.550, lambda=0.07781 
## + Fold2.Rep2: alpha=0.775, lambda=0.07781 
## - Fold2.Rep2: alpha=0.775, lambda=0.07781 
## + Fold2.Rep2: alpha=1.000, lambda=0.07781 
## - Fold2.Rep2: alpha=1.000, lambda=0.07781 
## + Fold3.Rep2: alpha=0.100, lambda=0.07781 
## - Fold3.Rep2: alpha=0.100, lambda=0.07781 
## + Fold3.Rep2: alpha=0.325, lambda=0.07781 
## - Fold3.Rep2: alpha=0.325, lambda=0.07781 
## + Fold3.Rep2: alpha=0.550, lambda=0.07781 
## - Fold3.Rep2: alpha=0.550, lambda=0.07781 
## + Fold3.Rep2: alpha=0.775, lambda=0.07781 
## - Fold3.Rep2: alpha=0.775, lambda=0.07781 
## + Fold3.Rep2: alpha=1.000, lambda=0.07781 
## - Fold3.Rep2: alpha=1.000, lambda=0.07781 
## + Fold1.Rep3: alpha=0.100, lambda=0.07781 
## - Fold1.Rep3: alpha=0.100, lambda=0.07781 
## + Fold1.Rep3: alpha=0.325, lambda=0.07781 
## - Fold1.Rep3: alpha=0.325, lambda=0.07781 
## + Fold1.Rep3: alpha=0.550, lambda=0.07781 
## - Fold1.Rep3: alpha=0.550, lambda=0.07781 
## + Fold1.Rep3: alpha=0.775, lambda=0.07781 
## - Fold1.Rep3: alpha=0.775, lambda=0.07781 
## + Fold1.Rep3: alpha=1.000, lambda=0.07781 
## - Fold1.Rep3: alpha=1.000, lambda=0.07781 
## + Fold2.Rep3: alpha=0.100, lambda=0.07781 
## - Fold2.Rep3: alpha=0.100, lambda=0.07781 
## + Fold2.Rep3: alpha=0.325, lambda=0.07781 
## - Fold2.Rep3: alpha=0.325, lambda=0.07781 
## + Fold2.Rep3: alpha=0.550, lambda=0.07781 
## - Fold2.Rep3: alpha=0.550, lambda=0.07781 
## + Fold2.Rep3: alpha=0.775, lambda=0.07781 
## - Fold2.Rep3: alpha=0.775, lambda=0.07781 
## + Fold2.Rep3: alpha=1.000, lambda=0.07781 
## - Fold2.Rep3: alpha=1.000, lambda=0.07781 
## + Fold3.Rep3: alpha=0.100, lambda=0.07781 
## - Fold3.Rep3: alpha=0.100, lambda=0.07781 
## + Fold3.Rep3: alpha=0.325, lambda=0.07781 
## - Fold3.Rep3: alpha=0.325, lambda=0.07781 
## + Fold3.Rep3: alpha=0.550, lambda=0.07781 
## - Fold3.Rep3: alpha=0.550, lambda=0.07781 
## + Fold3.Rep3: alpha=0.775, lambda=0.07781 
## - Fold3.Rep3: alpha=0.775, lambda=0.07781 
## + Fold3.Rep3: alpha=1.000, lambda=0.07781 
## - Fold3.Rep3: alpha=1.000, lambda=0.07781 
## Aggregating results
## Selecting tuning parameters
## Fitting alpha = 0.775, lambda = 0.00361 on full training set
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.data.training_0-2.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.data.training_0-3.png) 

```
##             Length Class      Mode     
## a0            89   -none-     numeric  
## beta        5073   dgCMatrix  S4       
## df            89   -none-     numeric  
## dim            2   -none-     numeric  
## lambda        89   -none-     numeric  
## dev.ratio     89   -none-     numeric  
## nulldev        1   -none-     numeric  
## npasses        1   -none-     numeric  
## jerr           1   -none-     numeric  
## offset         1   -none-     logical  
## classnames     2   -none-     character
## call           5   -none-     call     
## nobs           1   -none-     numeric  
## lambdaOpt      1   -none-     numeric  
## xNames        57   -none-     character
## problemType    1   -none-     character
## tuneValue      2   data.frame list     
## obsLevels      2   -none-     character
## [1] "min lambda > lambdaOpt:"
##                         (Intercept)                 NDSS.my.fctr#Mltmd# 
##                         -6.88661237                         -0.69017131 
##           NDSS.my.fctr#Opnn#RmFrDbt        NDSS.my.fctr#Opnn#ThPblcEdtr 
##                         -2.75091575                          3.16599950 
##             NDSS.my.fctr#U.S.#Edctn      NDSS.my.fctrBsnss#BsnssDy#Dlbk 
##                         -1.62723326                         -0.13226493 
## NDSS.my.fctrBsnss#BsnssDy#SmllBsnss     NDSS.my.fctrBsnss#Crsswrds/Gms# 
##                         -0.57337620                          3.30038315 
##          NDSS.my.fctrBsnss#Tchnlgy#              NDSS.my.fctrCltr#Arts# 
##                          0.49759684                         -0.02664973 
##              NDSS.my.fctrFrgn#Wrld#        NDSS.my.fctrFrgn#Wrld#AsPcfc 
##                         -0.84857323                         -1.38630180 
##              NDSS.my.fctrOpEd#Opnn#              NDSS.my.fctrScnc#Hlth# 
##                          3.78134805                          2.76279694 
##             NDSS.my.fctrStyls##Fshn             NDSS.my.fctrStyls#U.S.# 
##                         -1.04374095                          2.53368117 
##                 NDSS.my.fctrTStyl##              NDSS.my.fctrTrvl#Trvl# 
##                         -1.17568321                         -0.40646984 
##                  NDSS.my.fctrmyOthr            PubDate.date.fctr(13,19] 
##                         -0.92034972                         -0.03969260 
##          PubDate.day.minutes.poly.1          PubDate.day.minutes.poly.2 
##                         14.45900129                         15.50466097 
##          PubDate.day.minutes.poly.3          PubDate.day.minutes.poly.4 
##                          2.12223169                          4.01900421 
##        PubDate.hour.fctr(7.67,15.3]                PubDate.last16.log1p 
##                          0.13654098                          0.09222896 
##                PubDate.last32.log1p      PubDate.minute.fctr(29.5,44.2] 
##                          0.01367946                         -0.13143954 
##                PubDate.month.fctr11      PubDate.second.fctr(44.2,59.1] 
##                         -0.05213495                         -0.05593654 
##                 PubDate.wkday.fctr1                 PubDate.wkday.fctr5 
##                          0.08523292                         -0.11897940 
##                 PubDate.wkday.fctr6                       PubDate.wkend 
##                         -0.07242575                          0.30994932 
##                     WordCount.log1p                     WordCount.root2 
##                          0.35761159                          0.05338155 
## [1] "max lambda < lambdaOpt:"
##                         (Intercept)                 NDSS.my.fctr#Mltmd# 
##                        -7.024434156                        -0.773459716 
##           NDSS.my.fctr#Opnn#RmFrDbt        NDSS.my.fctr#Opnn#ThPblcEdtr 
##                        -2.859705647                         3.170283546 
##             NDSS.my.fctr#U.S.#Edctn      NDSS.my.fctrBsnss#BsnssDy#Dlbk 
##                        -1.742809921                        -0.166521901 
## NDSS.my.fctrBsnss#BsnssDy#SmllBsnss     NDSS.my.fctrBsnss#Crsswrds/Gms# 
##                        -0.624120498                         3.290504426 
##          NDSS.my.fctrBsnss#Tchnlgy#              NDSS.my.fctrCltr#Arts# 
##                         0.487995233                        -0.068894966 
##              NDSS.my.fctrFrgn#Wrld#        NDSS.my.fctrFrgn#Wrld#AsPcfc 
##                        -0.958508086                        -1.466854531 
##              NDSS.my.fctrOpEd#Opnn#              NDSS.my.fctrScnc#Hlth# 
##                         3.777480362                         2.750599151 
##             NDSS.my.fctrStyls##Fshn             NDSS.my.fctrStyls#U.S.# 
##                        -1.138065219                         2.520151641 
##                 NDSS.my.fctrTStyl##              NDSS.my.fctrTrvl#Trvl# 
##                        -1.236297003                        -0.481831991 
##                  NDSS.my.fctrmyOthr             PubDate.date.fctr(7,13] 
##                        -1.041571187                         0.002075209 
##            PubDate.date.fctr(13,19]          PubDate.day.minutes.poly.1 
##                        -0.048106489                        14.872950090 
##          PubDate.day.minutes.poly.2          PubDate.day.minutes.poly.3 
##                        16.952838409                         2.192634396 
##          PubDate.day.minutes.poly.4        PubDate.hour.fctr(7.67,15.3] 
##                         3.288655447                         0.178973137 
##                PubDate.last16.log1p                PubDate.last32.log1p 
##                         0.098282543                         0.015961358 
##      PubDate.minute.fctr(29.5,44.2]                PubDate.month.fctr11 
##                        -0.140081948                        -0.062964474 
##      PubDate.second.fctr(44.2,59.1]                 PubDate.wkday.fctr1 
##                        -0.064448103                         0.093086616 
##                 PubDate.wkday.fctr2                 PubDate.wkday.fctr5 
##                        -0.010913486                        -0.128522776 
##                 PubDate.wkday.fctr6                       PubDate.wkend 
##                        -0.110017248                         0.318726398 
##                     WordCount.log1p                     WordCount.root2 
##                         0.368532055                         0.053251505
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.data.training_0-4.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.data.training_0-5.png) 

```
##          Prediction
## Reference    N    Y
##         N 4981  458
##         Y  228  865
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   8.949786e-01   6.523492e-01   8.872900e-01   9.023124e-01   8.326699e-01 
## AccuracyPValue  McnemarPValue 
##   1.409434e-46   2.264774e-18 
##                  id
## 1 Final##rcv#glmnet
##                                                                                                                                                                                                                                                                                                                                                                                                                                                                feats
## 1 WordCount.root2,WordCount.log1p,NDSS.my.fctr,PubDate.day.minutes.poly.1,PubDate.hour.fctr,PubDate.wkend,PubDate.day.minutes.poly.4,PubDate.day.minutes.poly.2,PubDate.last4.log1p,PubDate.last2.log1p,PubDate.last8.log1p,PubDate.last16.log1p,PubDate.day.minutes.poly.3,PubDate.month.fctr,PubDate.juliandate,.rnorm,PubDate.last32.log1p,PubDate.date.fctr,PubDate.second.fctr,PubDate.minute.fctr,PubDate.wkday.fctr,WordCount.nexp,PubDate.day.minutes.poly.5
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1              25                     36.624                 0.817
##   max.AUCpROC.fit max.Sens.fit max.Spec.fit max.AUCROCR.fit
## 1        0.801372      0.96139    0.6413541       0.9355675
##   opt.prob.threshold.fit max.f.score.fit max.Accuracy.fit
## 1                    0.2       0.7160596        0.9065624
##   max.AccuracyLower.fit max.AccuracyUpper.fit max.Kappa.fit
## 1               0.88729             0.9023124     0.6401418
##   max.AccuracySD.fit max.KappaSD.fit
## 1        0.005702198      0.02639621
```

```r
rm(ret_lst)
glb_chunks_df <- myadd_chunk(glb_chunks_df, "fit.data.training", major.inc=FALSE)
```

```
##                label step_major step_minor label_minor     bgn     end
## 14 fit.data.training          7          0           0 192.948 258.677
## 15 fit.data.training          7          1           1 258.677      NA
##    elapsed
## 14  65.729
## 15      NA
```


```r
#stop(here"); glb2Sav()
if (glb_is_classification && glb_is_binomial) 
    prob_threshold <- glb_models_df[glb_models_df$id == glb_sel_mdl_id,
                                        "opt.prob.threshold.OOB"] else 
    prob_threshold <- NULL

if (grepl("Ensemble", glb_fin_mdl_id)) {
    # Get predictions for each model in ensemble; Outliers that have been moved to OOB might not have been predicted yet
    mdlEnsembleComps <- unlist(str_split(subset(glb_models_df, 
                                                id == glb_fin_mdl_id)$feats, ","))
    if (glb_is_classification && glb_is_binomial)
        mdlEnsembleComps <- gsub("\\.prob$", "", mdlEnsembleComps)
    mdlEnsembleComps <- gsub(paste0("^", 
                        gsub(".", "\\.", mygetPredictIds(glb_rsp_var)$value, fixed = TRUE)),
                             "", mdlEnsembleComps)
    for (mdl_id in mdlEnsembleComps) {
        glbObsTrn <- glb_get_predictions(df = glbObsTrn, mdl_id = mdl_id, 
                                            rsp_var = glb_rsp_var,
                                            prob_threshold_def = prob_threshold)
        glbObsNew <- glb_get_predictions(df = glbObsNew, mdl_id = mdl_id, 
                                            rsp_var = glb_rsp_var,
                                            prob_threshold_def = prob_threshold)
    }    
}
glbObsTrn <- glb_get_predictions(df = glbObsTrn, mdl_id = glb_fin_mdl_id, 
                                     rsp_var = glb_rsp_var,
                                    prob_threshold_def = prob_threshold)
```

```
## Warning in glb_get_predictions(df = glbObsTrn, mdl_id = glb_fin_mdl_id, :
## Using default probability threshold: 0.1
```

```r
glb_featsimp_df <- myget_feats_importance(mdl=glb_fin_mdl,
                                          featsimp_df=glb_featsimp_df)
glb_featsimp_df[, paste0(glb_fin_mdl_id, ".imp")] <- glb_featsimp_df$imp
print(glb_featsimp_df)
```

```
##                                     All.X##rcv#glmnet.imp        imp
## PubDate.day.minutes.poly.2                       5.307609 100.000000
## PubDate.day.minutes.poly.1                     100.000000  90.409031
## NDSS.my.fctrOpEd#Opnn#                          68.562523  33.933845
## PubDate.day.minutes.poly.4                       5.307609  32.183176
## NDSS.my.fctrBsnss#Crsswrds/Gms#                 63.559761  31.442333
## NDSS.my.fctr#Opnn#ThPblcEdtr                    48.417943  30.810953
## NDSS.my.fctrScnc#Hlth#                          54.307908  28.675494
## NDSS.my.fctrStyls#U.S.#                         50.421526  27.494881
## PubDate.day.minutes.poly.3                       5.307609  25.727586
## NDSS.my.fctrBsnss#Tchnlgy#                       5.690755  17.067558
## WordCount.log1p                                  6.701272  16.433441
## PubDate.wkend                                    5.345889  16.180212
## PubDate.hour.fctr(7.67,15.3]                     5.307609  15.428347
## PubDate.last16.log1p                             5.307609  15.052359
## PubDate.wkday.fctr1                              5.307609  15.023834
## WordCount.root2                                  6.021636  14.827828
## PubDate.last32.log1p                             5.307609  14.634050
## PubDate.date.fctr(7,13]                          5.307609  14.563041
## .rnorm                                           5.307609  14.554558
## NDSS.my.fctrCltr##                               5.307609  14.554558
## NDSS.my.fctrMtr#N.Y./Rgn#                        5.307609  14.554558
## PubDate.date.fctr(19,25]                         5.307609  14.554558
## PubDate.date.fctr(25,31]                         5.307609  14.554558
## PubDate.day.minutes.poly.5                       5.307609  14.554558
## PubDate.hour.fctr(15.3,23]                       5.307609  14.554558
## PubDate.juliandate                               5.307609  14.554558
## PubDate.last2.log1p                              5.307609  14.554558
## PubDate.last4.log1p                              5.307609  14.554558
## PubDate.last8.log1p                              5.307609  14.554558
## PubDate.minute.fctr(14.8,29.5]                   5.307609  14.554558
## PubDate.minute.fctr(44.2,59.1]                   5.307609  14.554558
## PubDate.month.fctr10                             5.307609  14.554558
## PubDate.month.fctr12                             5.307609  14.554558
## PubDate.second.fctr(14.8,29.5]                   5.307609  14.554558
## PubDate.second.fctr(29.5,44.2]                   5.307609  14.554558
## PubDate.wkday.fctr3                              5.307609  14.554558
## PubDate.wkday.fctr4                              5.307609  14.554558
## WordCount.nexp                                   5.307609  14.554558
## PubDate.wkday.fctr2                              5.307609  14.509947
## PubDate.date.fctr(13,19]                         5.307609  14.316575
## NDSS.my.fctrCltr#Arts#                           5.307609  14.245181
## PubDate.month.fctr11                             5.307609  14.242882
## PubDate.second.fctr(44.2,59.1]                   5.307609  14.232858
## PubDate.wkday.fctr6                              5.307609  14.029412
## PubDate.wkday.fctr5                              5.307609  13.905284
## PubDate.minute.fctr(29.5,44.2]                   5.307609  13.845057
## NDSS.my.fctrBsnss#BsnssDy#Dlbk                   5.307609  13.736119
## NDSS.my.fctrTrvl#Trvl#                           5.307609  12.161655
## NDSS.my.fctrBsnss#BsnssDy#SmllBsnss              5.307609  11.406198
## NDSS.my.fctr#Mltmd#                              5.307609  10.674108
## NDSS.my.fctrFrgn#Wrld#                           5.307609   9.752718
## NDSS.my.fctrmyOthr                               5.307609   9.338429
## NDSS.my.fctrStyls##Fshn                          5.307609   8.815485
## NDSS.my.fctrTStyl##                              2.967145   8.276532
## NDSS.my.fctrFrgn#Wrld#AsPcfc                     5.307609   7.114733
## NDSS.my.fctr#U.S.#Edctn                          5.307609   5.735793
## NDSS.my.fctr#Opnn#RmFrDbt                        0.000000   0.000000
##                                     Final##rcv#glmnet.imp
## PubDate.day.minutes.poly.2                     100.000000
## PubDate.day.minutes.poly.1                      90.409031
## NDSS.my.fctrOpEd#Opnn#                          33.933845
## PubDate.day.minutes.poly.4                      32.183176
## NDSS.my.fctrBsnss#Crsswrds/Gms#                 31.442333
## NDSS.my.fctr#Opnn#ThPblcEdtr                    30.810953
## NDSS.my.fctrScnc#Hlth#                          28.675494
## NDSS.my.fctrStyls#U.S.#                         27.494881
## PubDate.day.minutes.poly.3                      25.727586
## NDSS.my.fctrBsnss#Tchnlgy#                      17.067558
## WordCount.log1p                                 16.433441
## PubDate.wkend                                   16.180212
## PubDate.hour.fctr(7.67,15.3]                    15.428347
## PubDate.last16.log1p                            15.052359
## PubDate.wkday.fctr1                             15.023834
## WordCount.root2                                 14.827828
## PubDate.last32.log1p                            14.634050
## PubDate.date.fctr(7,13]                         14.563041
## .rnorm                                          14.554558
## NDSS.my.fctrCltr##                              14.554558
## NDSS.my.fctrMtr#N.Y./Rgn#                       14.554558
## PubDate.date.fctr(19,25]                        14.554558
## PubDate.date.fctr(25,31]                        14.554558
## PubDate.day.minutes.poly.5                      14.554558
## PubDate.hour.fctr(15.3,23]                      14.554558
## PubDate.juliandate                              14.554558
## PubDate.last2.log1p                             14.554558
## PubDate.last4.log1p                             14.554558
## PubDate.last8.log1p                             14.554558
## PubDate.minute.fctr(14.8,29.5]                  14.554558
## PubDate.minute.fctr(44.2,59.1]                  14.554558
## PubDate.month.fctr10                            14.554558
## PubDate.month.fctr12                            14.554558
## PubDate.second.fctr(14.8,29.5]                  14.554558
## PubDate.second.fctr(29.5,44.2]                  14.554558
## PubDate.wkday.fctr3                             14.554558
## PubDate.wkday.fctr4                             14.554558
## WordCount.nexp                                  14.554558
## PubDate.wkday.fctr2                             14.509947
## PubDate.date.fctr(13,19]                        14.316575
## NDSS.my.fctrCltr#Arts#                          14.245181
## PubDate.month.fctr11                            14.242882
## PubDate.second.fctr(44.2,59.1]                  14.232858
## PubDate.wkday.fctr6                             14.029412
## PubDate.wkday.fctr5                             13.905284
## PubDate.minute.fctr(29.5,44.2]                  13.845057
## NDSS.my.fctrBsnss#BsnssDy#Dlbk                  13.736119
## NDSS.my.fctrTrvl#Trvl#                          12.161655
## NDSS.my.fctrBsnss#BsnssDy#SmllBsnss             11.406198
## NDSS.my.fctr#Mltmd#                             10.674108
## NDSS.my.fctrFrgn#Wrld#                           9.752718
## NDSS.my.fctrmyOthr                               9.338429
## NDSS.my.fctrStyls##Fshn                          8.815485
## NDSS.my.fctrTStyl##                              8.276532
## NDSS.my.fctrFrgn#Wrld#AsPcfc                     7.114733
## NDSS.my.fctr#U.S.#Edctn                          5.735793
## NDSS.my.fctr#Opnn#RmFrDbt                        0.000000
```

```r
if (glb_is_classification && glb_is_binomial)
    glb_analytics_diag_plots(obs_df=glbObsTrn, mdl_id=glb_fin_mdl_id, 
            prob_threshold=glb_models_df[glb_models_df$id == glb_sel_mdl_id, 
                                         "opt.prob.threshold.OOB"]) else
    glb_analytics_diag_plots(obs_df=glbObsTrn, mdl_id=glb_fin_mdl_id)                  
```

```
## Warning in glb_analytics_diag_plots(obs_df = glbObsTrn, mdl_id =
## glb_fin_mdl_id, : Limiting important feature scatter plots to 5 out of 23
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.data.training_1-1.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.data.training_1-2.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.data.training_1-3.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.data.training_1-4.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.data.training_1-5.png) 

```
## [1] "Min/Max Boundaries: "
##   UniqueID Pplr.fctr Pplr.fctr.All.X..rcv.glmnet.prob
## 1     1065         N                               NA
## 2     4168         N                        0.0236808
## 3     5647         N                        0.1062812
## 4      302         N                               NA
##   Pplr.fctr.All.X..rcv.glmnet Pplr.fctr.All.X..rcv.glmnet.err
## 1                        <NA>                              NA
## 2                           N                           FALSE
## 3                           Y                            TRUE
## 4                        <NA>                              NA
##   Pplr.fctr.All.X..rcv.glmnet.err.abs Pplr.fctr.All.X..rcv.glmnet.is.acc
## 1                                  NA                                 NA
## 2                           0.0236808                               TRUE
## 3                           0.1062812                              FALSE
## 4                                  NA                                 NA
##   Pplr.fctr.Final..rcv.glmnet.prob Pplr.fctr.Final..rcv.glmnet
## 1                       0.02895198                           N
## 2                       0.00862190                           N
## 3                       0.11218487                           Y
## 4                       0.42886804                           Y
##   Pplr.fctr.Final..rcv.glmnet.err Pplr.fctr.Final..rcv.glmnet.err.abs
## 1                           FALSE                          0.02895198
## 2                           FALSE                          0.00862190
## 3                            TRUE                          0.11218487
## 4                            TRUE                          0.42886804
##   Pplr.fctr.Final..rcv.glmnet.is.acc Pplr.fctr.Final..rcv.glmnet.accurate
## 1                               TRUE                                 TRUE
## 2                               TRUE                                 TRUE
## 3                              FALSE                                FALSE
## 4                              FALSE                                FALSE
##   Pplr.fctr.Final..rcv.glmnet.error .label
## 1                        0.00000000   1065
## 2                        0.00000000   4168
## 3                        0.01218487   5647
## 4                        0.32886804    302
## [1] "Inaccurate: "
##   UniqueID Pplr.fctr Pplr.fctr.All.X..rcv.glmnet.prob
## 1     2182         Y                       0.02937511
## 2     4352         Y                               NA
## 3     4721         Y                       0.05500531
## 4     1696         Y                       0.02684498
## 5     5486         Y                               NA
## 6      364         Y                       0.03487556
##   Pplr.fctr.All.X..rcv.glmnet Pplr.fctr.All.X..rcv.glmnet.err
## 1                           N                            TRUE
## 2                        <NA>                              NA
## 3                           N                            TRUE
## 4                           N                            TRUE
## 5                        <NA>                              NA
## 6                           N                            TRUE
##   Pplr.fctr.All.X..rcv.glmnet.err.abs Pplr.fctr.All.X..rcv.glmnet.is.acc
## 1                           0.9706249                              FALSE
## 2                                  NA                                 NA
## 3                           0.9449947                              FALSE
## 4                           0.9731550                              FALSE
## 5                                  NA                                 NA
## 6                           0.9651244                              FALSE
##   Pplr.fctr.Final..rcv.glmnet.prob Pplr.fctr.Final..rcv.glmnet
## 1                      0.007066632                           N
## 2                      0.011579123                           N
## 3                      0.014803392                           N
## 4                      0.014856448                           N
## 5                      0.017198535                           N
## 6                      0.019459474                           N
##   Pplr.fctr.Final..rcv.glmnet.err Pplr.fctr.Final..rcv.glmnet.err.abs
## 1                            TRUE                           0.9929334
## 2                            TRUE                           0.9884209
## 3                            TRUE                           0.9851966
## 4                            TRUE                           0.9851436
## 5                            TRUE                           0.9828015
## 6                            TRUE                           0.9805405
##   Pplr.fctr.Final..rcv.glmnet.is.acc Pplr.fctr.Final..rcv.glmnet.accurate
## 1                              FALSE                                FALSE
## 2                              FALSE                                FALSE
## 3                              FALSE                                FALSE
## 4                              FALSE                                FALSE
## 5                              FALSE                                FALSE
## 6                              FALSE                                FALSE
##   Pplr.fctr.Final..rcv.glmnet.error
## 1                       -0.09293337
## 2                       -0.08842088
## 3                       -0.08519661
## 4                       -0.08514355
## 5                       -0.08280147
## 6                       -0.08054053
##      UniqueID Pplr.fctr Pplr.fctr.All.X..rcv.glmnet.prob
## 191       988         N                       0.07208817
## 221      3609         N                       0.11095536
## 382      4016         N                       0.06234938
## 413      3378         N                               NA
## 909      1805         N                       0.51677328
## 1110     6511         N                               NA
##      Pplr.fctr.All.X..rcv.glmnet Pplr.fctr.All.X..rcv.glmnet.err
## 191                            N                           FALSE
## 221                            Y                            TRUE
## 382                            N                           FALSE
## 413                         <NA>                              NA
## 909                            Y                            TRUE
## 1110                        <NA>                              NA
##      Pplr.fctr.All.X..rcv.glmnet.err.abs
## 191                           0.07208817
## 221                           0.11095536
## 382                           0.06234938
## 413                                   NA
## 909                           0.51677328
## 1110                                  NA
##      Pplr.fctr.All.X..rcv.glmnet.is.acc Pplr.fctr.Final..rcv.glmnet.prob
## 191                                TRUE                        0.1117584
## 221                               FALSE                        0.1158289
## 382                                TRUE                        0.1346807
## 413                                  NA                        0.1377206
## 909                               FALSE                        0.3730036
## 1110                                 NA                        0.7233569
##      Pplr.fctr.Final..rcv.glmnet Pplr.fctr.Final..rcv.glmnet.err
## 191                            Y                            TRUE
## 221                            Y                            TRUE
## 382                            Y                            TRUE
## 413                            Y                            TRUE
## 909                            Y                            TRUE
## 1110                           Y                            TRUE
##      Pplr.fctr.Final..rcv.glmnet.err.abs
## 191                            0.1117584
## 221                            0.1158289
## 382                            0.1346807
## 413                            0.1377206
## 909                            0.3730036
## 1110                           0.7233569
##      Pplr.fctr.Final..rcv.glmnet.is.acc
## 191                               FALSE
## 221                               FALSE
## 382                               FALSE
## 413                               FALSE
## 909                               FALSE
## 1110                              FALSE
##      Pplr.fctr.Final..rcv.glmnet.accurate
## 191                                 FALSE
## 221                                 FALSE
## 382                                 FALSE
## 413                                 FALSE
## 909                                 FALSE
## 1110                                FALSE
##      Pplr.fctr.Final..rcv.glmnet.error
## 191                         0.01175842
## 221                         0.01582887
## 382                         0.03468074
## 413                         0.03772058
## 909                         0.27300359
## 1110                        0.62335692
##      UniqueID Pplr.fctr Pplr.fctr.All.X..rcv.glmnet.prob
## 1188      770         N                               NA
## 1189     2179         N                        0.9167004
## 1190      472         N                               NA
## 1191     2995         N                               NA
## 1192     1612         N                        0.9192335
## 1193     1448         N                               NA
##      Pplr.fctr.All.X..rcv.glmnet Pplr.fctr.All.X..rcv.glmnet.err
## 1188                        <NA>                              NA
## 1189                           Y                            TRUE
## 1190                        <NA>                              NA
## 1191                        <NA>                              NA
## 1192                           Y                            TRUE
## 1193                        <NA>                              NA
##      Pplr.fctr.All.X..rcv.glmnet.err.abs
## 1188                                  NA
## 1189                           0.9167004
## 1190                                  NA
## 1191                                  NA
## 1192                           0.9192335
## 1193                                  NA
##      Pplr.fctr.All.X..rcv.glmnet.is.acc Pplr.fctr.Final..rcv.glmnet.prob
## 1188                                 NA                        0.9488795
## 1189                              FALSE                        0.9498766
## 1190                                 NA                        0.9568880
## 1191                                 NA                        0.9602272
## 1192                              FALSE                        0.9635908
## 1193                                 NA                        0.9643650
##      Pplr.fctr.Final..rcv.glmnet Pplr.fctr.Final..rcv.glmnet.err
## 1188                           Y                            TRUE
## 1189                           Y                            TRUE
## 1190                           Y                            TRUE
## 1191                           Y                            TRUE
## 1192                           Y                            TRUE
## 1193                           Y                            TRUE
##      Pplr.fctr.Final..rcv.glmnet.err.abs
## 1188                           0.9488795
## 1189                           0.9498766
## 1190                           0.9568880
## 1191                           0.9602272
## 1192                           0.9635908
## 1193                           0.9643650
##      Pplr.fctr.Final..rcv.glmnet.is.acc
## 1188                              FALSE
## 1189                              FALSE
## 1190                              FALSE
## 1191                              FALSE
## 1192                              FALSE
## 1193                              FALSE
##      Pplr.fctr.Final..rcv.glmnet.accurate
## 1188                                FALSE
## 1189                                FALSE
## 1190                                FALSE
## 1191                                FALSE
## 1192                                FALSE
## 1193                                FALSE
##      Pplr.fctr.Final..rcv.glmnet.error
## 1188                         0.8488795
## 1189                         0.8498766
## 1190                         0.8568880
## 1191                         0.8602272
## 1192                         0.8635908
## 1193                         0.8643650
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.data.training_1-6.png) 

```r
dsp_feats_vctr <- c(NULL)
for(var in grep(".imp", names(glb_feats_df), fixed=TRUE, value=TRUE))
    dsp_feats_vctr <- union(dsp_feats_vctr, 
                            glb_feats_df[!is.na(glb_feats_df[, var]), "id"])

# print(glbObsTrn[glbObsTrn$UniqueID %in% FN_OOB_ids, 
#                     grep(glb_rsp_var, names(glbObsTrn), value=TRUE)])

print(setdiff(names(glbObsTrn), names(glbObsAll)))
```

```
## [1] "Pplr.fctr.Final..rcv.glmnet.prob"   
## [2] "Pplr.fctr.Final..rcv.glmnet"        
## [3] "Pplr.fctr.Final..rcv.glmnet.err"    
## [4] "Pplr.fctr.Final..rcv.glmnet.err.abs"
## [5] "Pplr.fctr.Final..rcv.glmnet.is.acc"
```

```r
for (col in setdiff(names(glbObsTrn), names(glbObsAll)))
    # Merge or cbind ?
    glbObsAll[glbObsAll$.src == "Train", col] <- glbObsTrn[, col]

print(setdiff(names(glbObsFit), names(glbObsAll)))
```

```
## character(0)
```

```r
print(setdiff(names(glbObsOOB), names(glbObsAll)))
```

```
## character(0)
```

```r
for (col in setdiff(names(glbObsOOB), names(glbObsAll)))
    # Merge or cbind ?
    glbObsAll[glbObsAll$.lcn == "OOB", col] <- glbObsOOB[, col]
    
print(setdiff(names(glbObsNew), names(glbObsAll)))
```

```
## character(0)
```

```r
if (glb_save_envir)
    save(glb_feats_df, glbObsAll, 
         #glbObsTrn, glbObsFit, glbObsOOB, glbObsNew,
         glb_models_df, dsp_models_df, glb_models_lst, glb_model_type,
         glb_sel_mdl, glb_sel_mdl_id,
         glb_fin_mdl, glb_fin_mdl_id,
        file=paste0(glb_out_pfx, "dsk.RData"))

replay.petrisim(pn=glb_analytics_pn, 
    replay.trans=(glb_analytics_avl_objs <- c(glb_analytics_avl_objs, 
        "data.training.all.prediction","model.final")), flip_coord=TRUE)
```

```
## time	trans	 "bgn " "fit.data.training.all " "predict.data.new " "end " 
## 0.0000 	multiple enabled transitions:  data.training.all data.new model.selected 	firing:  data.training.all 
## 1.0000 	 1 	 2 1 0 0 
## 1.0000 	multiple enabled transitions:  data.training.all data.new model.selected model.final data.training.all.prediction 	firing:  data.new 
## 2.0000 	 2 	 1 1 1 0 
## 2.0000 	multiple enabled transitions:  data.training.all data.new model.selected model.final data.training.all.prediction data.new.prediction 	firing:  model.selected 
## 3.0000 	 3 	 0 2 1 0 
## 3.0000 	multiple enabled transitions:  model.final data.training.all.prediction data.new.prediction 	firing:  data.training.all.prediction 
## 4.0000 	 5 	 0 1 1 1 
## 4.0000 	multiple enabled transitions:  model.final data.training.all.prediction data.new.prediction 	firing:  model.final 
## 5.0000 	 4 	 0 0 2 1
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.data.training_1-7.png) 

```r
glb_chunks_df <- myadd_chunk(glb_chunks_df, "predict.data.new", major.inc=TRUE)
```

```
##                label step_major step_minor label_minor     bgn     end
## 15 fit.data.training          7          1           1 258.677 269.153
## 16  predict.data.new          8          0           0 269.154      NA
##    elapsed
## 15  10.476
## 16      NA
```

## Step `8.0: predict data new`

```
## Warning in glb_get_predictions(obs_df, mdl_id = glb_fin_mdl_id, rsp_var =
## glb_rsp_var, : Using default probability threshold: 0.1
```

```
## Warning in glb_get_predictions(obs_df, mdl_id = glb_fin_mdl_id, rsp_var =
## glb_rsp_var, : Using default probability threshold: 0.1
```

```
## Warning in zoo(rval[i], index(x)[i]): some methods for "zoo" objects do not
## work if the index entries in 'order.by' are not unique
```

```
## Warning in glb_analytics_diag_plots(obs_df = glbObsNew, mdl_id =
## glb_fin_mdl_id, : Limiting important feature scatter plots to 5 out of 23
```

```
## Warning: Removed 1870 rows containing missing values (geom_point).
```

```
## Warning: Removed 1870 rows containing missing values (geom_point).
```

![](NYTBlogs3_feat_PubDate_files/figure-html/predict.data.new-1.png) 

```
## Warning: Removed 1870 rows containing missing values (geom_point).
```

```
## Warning: Removed 1870 rows containing missing values (geom_point).
```

![](NYTBlogs3_feat_PubDate_files/figure-html/predict.data.new-2.png) 

```
## Warning: Removed 1870 rows containing missing values (geom_point).
```

```
## Warning: Removed 1870 rows containing missing values (geom_point).
```

![](NYTBlogs3_feat_PubDate_files/figure-html/predict.data.new-3.png) 

```
## Warning: Removed 1870 rows containing missing values (geom_point).
```

```
## Warning: Removed 1870 rows containing missing values (geom_point).
```

![](NYTBlogs3_feat_PubDate_files/figure-html/predict.data.new-4.png) 

```
## Warning: Removed 1870 rows containing missing values (geom_point).
```

```
## Warning: Removed 1870 rows containing missing values (geom_point).
```

![](NYTBlogs3_feat_PubDate_files/figure-html/predict.data.new-5.png) 

```
## NULL
```

```
## Loading required package: stringr
```

```
## [1] "ObsNew Prediction errors in categories:"
##    NDSS.my.fctr .n.Trn.N .n.Trn.Y .n.New.N .n.New.Y
## 5   #U.S.#Edctn      325       NA       87        2
## 10       Cltr##        1       NA       48       22
## .n.Trn.N .n.Trn.Y .n.New.N .n.New.Y 
##      326        0      135       24
```

```
## Loading required package: tidyr
## 
## Attaching package: 'tidyr'
## 
## The following object is masked from 'package:Matrix':
## 
##     expand
```

```
## [1] "OOBobs Pplr.fctr.All.X..rcv.glmnet Y: min < min of Train range: 2"
##      UniqueID Pplr.fctr.All.X..rcv.glmnet PubDate.day.minutes.poly.2
## 6435     6435                           Y                0.002105741
## 4223     4223                           Y               -0.008758791
##      WordCount.log1p WordCount.root2
## 6435         0.00000         0.00000
## 4223         7.27448        37.97368
##                                                    id      cor.y
## PubDate.day.minutes.poly.2 PubDate.day.minutes.poly.2 0.07097772
## WordCount.log1p                       WordCount.log1p 0.25431963
## WordCount.root2                       WordCount.root2 0.29212068
##                            exclude.as.feat  cor.y.abs      cor.high.X
## PubDate.day.minutes.poly.2           FALSE 0.07097772            <NA>
## WordCount.log1p                      FALSE 0.25431963 WordCount.root2
## WordCount.root2                      FALSE 0.29212068            <NA>
##                            freqRatio percentUnique zeroVar   nzv
## PubDate.day.minutes.poly.2  1.225490      18.08022   FALSE FALSE
## WordCount.log1p             2.315789      24.15799   FALSE FALSE
## WordCount.root2             2.315789      24.15799   FALSE FALSE
##                            is.cor.y.abs.low interaction.feat
## PubDate.day.minutes.poly.2            FALSE               NA
## WordCount.log1p                       FALSE               NA
## WordCount.root2                       FALSE               NA
##                            shapiro.test.p.value rsp_var_raw id_var rsp_var
## PubDate.day.minutes.poly.2         8.020999e-64       FALSE     NA      NA
## WordCount.log1p                    1.576866e-49       FALSE     NA      NA
## WordCount.root2                    4.556481e-30       FALSE     NA      NA
##                                     max          min max.Pplr.fctr.N
## PubDate.day.minutes.poly.2   0.04268445 -0.008758791      0.04268445
## WordCount.log1p              9.29771002  0.000000000      8.81966535
## WordCount.root2            104.46051886  0.000000000     82.24962006
##                            max.Pplr.fctr.Y min.Pplr.fctr.N min.Pplr.fctr.Y
## PubDate.day.minutes.poly.2      0.04254377    -0.008758791    -0.008758717
## WordCount.log1p                 9.29771002     0.000000000     1.945910149
## WordCount.root2               104.46051886     0.000000000     2.449489743
##                            max.Pplr.fctr.All.X..rcv.glmnet.N
## PubDate.day.minutes.poly.2                        0.04268445
## WordCount.log1p                                   7.06902343
## WordCount.root2                                  34.26368340
##                            max.Pplr.fctr.All.X..rcv.glmnet.Y
## PubDate.day.minutes.poly.2                        0.04254377
## WordCount.log1p                                   9.14088311
## WordCount.root2                                  96.58157174
##                            min.Pplr.fctr.All.X..rcv.glmnet.N
## PubDate.day.minutes.poly.2                      -0.008758672
## WordCount.log1p                                  0.000000000
## WordCount.root2                                  0.000000000
##                            min.Pplr.fctr.All.X..rcv.glmnet.Y
## PubDate.day.minutes.poly.2                      -0.008758791
## WordCount.log1p                                  0.000000000
## WordCount.root2                                  0.000000000
##                            max.Pplr.fctr.Final..rcv.glmnet.N
## PubDate.day.minutes.poly.2                        0.04268445
## WordCount.log1p                                   7.94093976
## WordCount.root2                                  53.00000000
##                            max.Pplr.fctr.Final..rcv.glmnet.Y
## PubDate.day.minutes.poly.2                        0.04254377
## WordCount.log1p                                   8.69232228
## WordCount.root2                                  77.17512553
##                            min.Pplr.fctr.Final..rcv.glmnet.N
## PubDate.day.minutes.poly.2                      -0.008758791
## WordCount.log1p                                  0.000000000
## WordCount.root2                                  0.000000000
##                            min.Pplr.fctr.Final..rcv.glmnet.Y
## PubDate.day.minutes.poly.2                      -0.008758672
## WordCount.log1p                                  1.609437912
## WordCount.root2                                  2.000000000
## [1] "OOBobs Pplr.fctr.All.X..rcv.glmnet Y: max > max of Train range: 5"
##      UniqueID Pplr.fctr.All.X..rcv.glmnet PubDate.day.minutes.poly.1
## 6528     6528                           Y                0.001809613
## 302       302                           Y                0.024722851
## 6435     6435                           Y               -0.013151170
## 3205     3205                           Y                0.002898990
## 6521     6521                           Y                0.013901702
##      PubDate.day.minutes.poly.3 PubDate.day.minutes.poly.4
## 6528               -0.004472838                0.007219141
## 302                 0.051829879                0.066100937
## 6435                0.010843928               -0.010140008
## 3205               -0.005881110                0.005605493
## 6521               -0.004412350               -0.013803062
##      PubDate.day.minutes.poly.5 PubDate.last16.log1p PubDate.last32.log1p
## 6528                0.006658270             11.95698             12.10953
## 302                 0.083442278             10.07706             10.33290
## 6435               -0.001891593             10.33864             10.86697
## 3205                0.008323368             11.85059             12.00544
## 6521               -0.012191759             11.68539             12.19622
##      PubDate.last8.log1p WordCount.nexp
## 6528           11.425547   0.000000e+00
## 302             9.741557   0.000000e+00
## 6435            8.695172   1.000000e+00
## 3205           11.443361   1.026188e-10
## 6521           10.414633   0.000000e+00
##                                                    id        cor.y
## PubDate.day.minutes.poly.1 PubDate.day.minutes.poly.1  0.156753478
## PubDate.day.minutes.poly.3 PubDate.day.minutes.poly.3  0.027983551
## PubDate.day.minutes.poly.4 PubDate.day.minutes.poly.4  0.073941394
## PubDate.day.minutes.poly.5 PubDate.day.minutes.poly.5 -0.055929231
## PubDate.last16.log1p             PubDate.last16.log1p  0.040735543
## PubDate.last32.log1p             PubDate.last32.log1p  0.003558081
## PubDate.last8.log1p               PubDate.last8.log1p  0.054458821
## WordCount.nexp                         WordCount.nexp -0.053208396
##                            exclude.as.feat   cor.y.abs          cor.high.X
## PubDate.day.minutes.poly.1           FALSE 0.156753478                <NA>
## PubDate.day.minutes.poly.3           FALSE 0.027983551                <NA>
## PubDate.day.minutes.poly.4           FALSE 0.073941394                <NA>
## PubDate.day.minutes.poly.5           FALSE 0.055929231                <NA>
## PubDate.last16.log1p                 FALSE 0.040735543                <NA>
## PubDate.last32.log1p                 FALSE 0.003558081                <NA>
## PubDate.last8.log1p                  FALSE 0.054458821 PubDate.last4.log1p
## WordCount.nexp                       FALSE 0.053208396                <NA>
##                            freqRatio percentUnique zeroVar   nzv
## PubDate.day.minutes.poly.1  1.225490      18.08022   FALSE FALSE
## PubDate.day.minutes.poly.3  1.225490      18.08022   FALSE FALSE
## PubDate.day.minutes.poly.4  1.225490      18.08022   FALSE FALSE
## PubDate.day.minutes.poly.5  1.225490      18.08022   FALSE FALSE
## PubDate.last16.log1p        3.200000      84.44581   FALSE FALSE
## PubDate.last32.log1p        8.000000      90.99816   FALSE FALSE
## PubDate.last8.log1p         1.142857      75.12247   FALSE FALSE
## WordCount.nexp             17.761364      11.32884   FALSE FALSE
##                            is.cor.y.abs.low interaction.feat
## PubDate.day.minutes.poly.1            FALSE               NA
## PubDate.day.minutes.poly.3            FALSE               NA
## PubDate.day.minutes.poly.4            FALSE               NA
## PubDate.day.minutes.poly.5            FALSE               NA
## PubDate.last16.log1p                  FALSE               NA
## PubDate.last32.log1p                   TRUE               NA
## PubDate.last8.log1p                   FALSE               NA
## WordCount.nexp                        FALSE               NA
##                            shapiro.test.p.value rsp_var_raw id_var rsp_var
## PubDate.day.minutes.poly.1         1.590362e-18       FALSE     NA      NA
## PubDate.day.minutes.poly.3         9.822405e-52       FALSE     NA      NA
## PubDate.day.minutes.poly.4         1.523136e-47       FALSE     NA      NA
## PubDate.day.minutes.poly.5         1.157500e-41       FALSE     NA      NA
## PubDate.last16.log1p               7.310334e-68       FALSE     NA      NA
## PubDate.last32.log1p               2.783236e-77       FALSE     NA      NA
## PubDate.last8.log1p                3.859176e-56       FALSE     NA      NA
## WordCount.nexp                     9.108805e-94       FALSE     NA      NA
##                                    max         min max.Pplr.fctr.N
## PubDate.day.minutes.poly.1  0.02475916 -0.02749464      0.02468654
## PubDate.day.minutes.poly.3  0.05215301 -0.04512497      0.05150779
## PubDate.day.minutes.poly.4  0.06677441 -0.01832740      0.06543120
## PubDate.day.minutes.poly.5  0.08471756 -0.02450918      0.08217780
## PubDate.last16.log1p       11.95698288  0.00000000     11.94531808
## PubDate.last32.log1p       12.32340669  0.00000000     12.21244232
## PubDate.last8.log1p        11.62246125  0.00000000     11.43577441
## WordCount.nexp              1.00000000  0.00000000      1.00000000
##                            max.Pplr.fctr.Y min.Pplr.fctr.N min.Pplr.fctr.Y
## PubDate.day.minutes.poly.1     0.024468663     -0.02749464     -0.02745833
## PubDate.day.minutes.poly.3     0.049597025     -0.04512497     -0.04482024
## PubDate.day.minutes.poly.4     0.061490534     -0.01832740     -0.01821959
## PubDate.day.minutes.poly.5     0.074814724     -0.02450918     -0.02362780
## PubDate.last16.log1p          11.877603300      0.00000000      0.00000000
## PubDate.last32.log1p          12.178408497      0.00000000      0.00000000
## PubDate.last8.log1p           11.394288315      0.00000000      0.00000000
## WordCount.nexp                 0.002478752      0.00000000      0.00000000
##                            max.Pplr.fctr.All.X..rcv.glmnet.N
## PubDate.day.minutes.poly.1                        0.02450498
## PubDate.day.minutes.poly.3                        0.04991290
## PubDate.day.minutes.poly.4                        0.06213811
## PubDate.day.minutes.poly.5                        0.07601554
## PubDate.last16.log1p                             11.86531323
## PubDate.last32.log1p                             12.21791228
## PubDate.last8.log1p                              11.62246125
## WordCount.nexp                                    1.00000000
##                            max.Pplr.fctr.All.X..rcv.glmnet.Y
## PubDate.day.minutes.poly.1                        0.02472285
## PubDate.day.minutes.poly.3                        0.05182988
## PubDate.day.minutes.poly.4                        0.06610094
## PubDate.day.minutes.poly.5                        0.08344228
## PubDate.last16.log1p                             11.95698288
## PubDate.last32.log1p                             12.19622431
## PubDate.last8.log1p                              11.44336100
## WordCount.nexp                                    1.00000000
##                            min.Pplr.fctr.All.X..rcv.glmnet.N
## PubDate.day.minutes.poly.1                       -0.02749464
## PubDate.day.minutes.poly.3                       -0.04512497
## PubDate.day.minutes.poly.4                       -0.01832685
## PubDate.day.minutes.poly.5                       -0.02450918
## PubDate.last16.log1p                              0.00000000
## PubDate.last32.log1p                              0.00000000
## PubDate.last8.log1p                               0.00000000
## WordCount.nexp                                    0.00000000
##                            min.Pplr.fctr.All.X..rcv.glmnet.Y
## PubDate.day.minutes.poly.1                       -0.02745833
## PubDate.day.minutes.poly.3                       -0.04482024
## PubDate.day.minutes.poly.4                       -0.01791547
## PubDate.day.minutes.poly.5                       -0.02362780
## PubDate.last16.log1p                              0.00000000
## PubDate.last32.log1p                              0.00000000
## PubDate.last8.log1p                               7.12447826
## WordCount.nexp                                    0.00000000
##                            max.Pplr.fctr.Final..rcv.glmnet.N
## PubDate.day.minutes.poly.1                        0.02432341
## PubDate.day.minutes.poly.3                        0.04834381
## PubDate.day.minutes.poly.4                        0.05893666
## PubDate.day.minutes.poly.5                        0.07011504
## PubDate.last16.log1p                             11.88113167
## PubDate.last32.log1p                             12.32340669
## PubDate.last8.log1p                              11.27955479
## WordCount.nexp                                    1.00000000
##                            max.Pplr.fctr.Final..rcv.glmnet.Y
## PubDate.day.minutes.poly.1                        0.02475916
## PubDate.day.minutes.poly.3                        0.05215301
## PubDate.day.minutes.poly.4                        0.06677441
## PubDate.day.minutes.poly.5                        0.08471756
## PubDate.last16.log1p                             11.84328621
## PubDate.last32.log1p                             12.30546086
## PubDate.last8.log1p                              11.33227851
## WordCount.nexp                                    0.01831564
##                            min.Pplr.fctr.Final..rcv.glmnet.N
## PubDate.day.minutes.poly.1                       -0.02749464
## PubDate.day.minutes.poly.3                       -0.04512497
## PubDate.day.minutes.poly.4                       -0.01832268
## PubDate.day.minutes.poly.5                       -0.02450918
## PubDate.last16.log1p                              8.10167775
## PubDate.last32.log1p                              8.86290830
## PubDate.last8.log1p                               7.06133437
## WordCount.nexp                                    0.00000000
##                            min.Pplr.fctr.Final..rcv.glmnet.Y
## PubDate.day.minutes.poly.1                       -0.02745833
## PubDate.day.minutes.poly.3                       -0.04482024
## PubDate.day.minutes.poly.4                       -0.01820339
## PubDate.day.minutes.poly.5                       -0.02362780
## PubDate.last16.log1p                              8.09223941
## PubDate.last32.log1p                              8.83579237
## PubDate.last8.log1p                               6.89162590
## WordCount.nexp                                    0.00000000
## [1] "OOBobs Pplr.fctr.All.X..rcv.glmnet N: max > max of Train range: 2"
##      UniqueID Pplr.fctr.All.X..rcv.glmnet PubDate.last32.log1p
## 5233     5233                           N             11.81975
## 6517     6517                           N             12.21791
##      PubDate.last8.log1p
## 5233           11.622461
## 6517            9.753188
##                                        id       cor.y exclude.as.feat
## PubDate.last32.log1p PubDate.last32.log1p 0.003558081           FALSE
## PubDate.last8.log1p   PubDate.last8.log1p 0.054458821           FALSE
##                        cor.y.abs          cor.high.X freqRatio
## PubDate.last32.log1p 0.003558081                <NA>  8.000000
## PubDate.last8.log1p  0.054458821 PubDate.last4.log1p  1.142857
##                      percentUnique zeroVar   nzv is.cor.y.abs.low
## PubDate.last32.log1p      90.99816   FALSE FALSE             TRUE
## PubDate.last8.log1p       75.12247   FALSE FALSE            FALSE
##                      interaction.feat shapiro.test.p.value rsp_var_raw
## PubDate.last32.log1p               NA         2.783236e-77       FALSE
## PubDate.last8.log1p                NA         3.859176e-56       FALSE
##                      id_var rsp_var      max min max.Pplr.fctr.N
## PubDate.last32.log1p     NA      NA 12.32341   0        12.21244
## PubDate.last8.log1p      NA      NA 11.62246   0        11.43577
##                      max.Pplr.fctr.Y min.Pplr.fctr.N min.Pplr.fctr.Y
## PubDate.last32.log1p        12.17841               0               0
## PubDate.last8.log1p         11.39429               0               0
##                      max.Pplr.fctr.All.X..rcv.glmnet.N
## PubDate.last32.log1p                          12.21791
## PubDate.last8.log1p                           11.62246
##                      max.Pplr.fctr.All.X..rcv.glmnet.Y
## PubDate.last32.log1p                          12.19622
## PubDate.last8.log1p                           11.44336
##                      min.Pplr.fctr.All.X..rcv.glmnet.N
## PubDate.last32.log1p                                 0
## PubDate.last8.log1p                                  0
##                      min.Pplr.fctr.All.X..rcv.glmnet.Y
## PubDate.last32.log1p                          0.000000
## PubDate.last8.log1p                           7.124478
##                      max.Pplr.fctr.Final..rcv.glmnet.N
## PubDate.last32.log1p                          12.32341
## PubDate.last8.log1p                           11.27955
##                      max.Pplr.fctr.Final..rcv.glmnet.Y
## PubDate.last32.log1p                          12.30546
## PubDate.last8.log1p                           11.33228
##                      min.Pplr.fctr.Final..rcv.glmnet.N
## PubDate.last32.log1p                          8.862908
## PubDate.last8.log1p                           7.061334
##                      min.Pplr.fctr.Final..rcv.glmnet.Y
## PubDate.last32.log1p                          8.835792
## PubDate.last8.log1p                           6.891626
## [1] "OOBobs total range outliers: 8"
## [1] "newobs Pplr.fctr.Final..rcv.glmnet N: max > max of Train range: 1186"
##      UniqueID Pplr.fctr.Final..rcv.glmnet PubDate.juliandate
## 6533     6533                           N                335
## 6540     6540                           N                335
## 6541     6541                           N                335
## 6542     6542                           N                335
## 6543     6543                           N                335
## 6545     6545                           N                335
##      PubDate.last32.log1p
## 6533            10.134321
## 6540             9.418817
## 6541             9.442800
## 6542             9.391411
## 6543             9.400878
## 6545             9.402860
##      UniqueID Pplr.fctr.Final..rcv.glmnet PubDate.juliandate
## 6629     6629                           N                335
## 7478     7478                           N                346
## 7575     7575                           N                349
## 7720     7720                           N                351
## 7907     7907                           N                352
## 8236     8236                           N                360
##      PubDate.last32.log1p
## 6629            12.219719
## 7478            10.902777
## 7575            10.133924
## 7720             9.276783
## 7907            10.825760
## 8236            11.729246
##      UniqueID Pplr.fctr.Final..rcv.glmnet PubDate.juliandate
## 8394     8394                           N                365
## 8395     8395                           N                365
## 8396     8396                           N                365
## 8398     8398                           N                365
## 8400     8400                           N                365
## 8401     8401                           N                365
##      PubDate.last32.log1p
## 8394             11.14224
## 8395             11.14428
## 8396             11.11651
## 8398             11.09410
## 8400             11.07855
## 8401             11.05398
##                                        id       cor.y exclude.as.feat
## PubDate.juliandate     PubDate.juliandate 0.014361075           FALSE
## PubDate.last32.log1p PubDate.last32.log1p 0.003558081           FALSE
##                        cor.y.abs         cor.high.X freqRatio
## PubDate.juliandate   0.014361075 PubDate.month.fctr   1.03252
## PubDate.last32.log1p 0.003558081               <NA>   8.00000
##                      percentUnique zeroVar   nzv is.cor.y.abs.low
## PubDate.juliandate        1.393141   FALSE FALSE            FALSE
## PubDate.last32.log1p     90.998163   FALSE FALSE             TRUE
##                      interaction.feat shapiro.test.p.value rsp_var_raw
## PubDate.juliandate                 NA         1.389406e-35       FALSE
## PubDate.last32.log1p               NA         2.783236e-77       FALSE
##                      id_var rsp_var       max min max.Pplr.fctr.N
## PubDate.juliandate       NA      NA 365.00000 244       334.00000
## PubDate.last32.log1p     NA      NA  12.32341   0        12.21244
##                      max.Pplr.fctr.Y min.Pplr.fctr.N min.Pplr.fctr.Y
## PubDate.juliandate         334.00000             244             244
## PubDate.last32.log1p        12.21791               0               0
##                      max.Pplr.fctr.All.X..rcv.glmnet.N
## PubDate.juliandate                           334.00000
## PubDate.last32.log1p                          12.21791
##                      max.Pplr.fctr.All.X..rcv.glmnet.Y
## PubDate.juliandate                           334.00000
## PubDate.last32.log1p                          12.19622
##                      min.Pplr.fctr.All.X..rcv.glmnet.N
## PubDate.juliandate                                 244
## PubDate.last32.log1p                                 0
##                      min.Pplr.fctr.All.X..rcv.glmnet.Y
## PubDate.juliandate                                 244
## PubDate.last32.log1p                                 0
##                      max.Pplr.fctr.Final..rcv.glmnet.N
## PubDate.juliandate                           365.00000
## PubDate.last32.log1p                          12.32341
##                      max.Pplr.fctr.Final..rcv.glmnet.Y
## PubDate.juliandate                           365.00000
## PubDate.last32.log1p                          12.30546
##                      min.Pplr.fctr.Final..rcv.glmnet.N
## PubDate.juliandate                          335.000000
## PubDate.last32.log1p                          8.862908
##                      min.Pplr.fctr.Final..rcv.glmnet.Y
## PubDate.juliandate                          335.000000
## PubDate.last32.log1p                          8.835792
## [1] "newobs Pplr.fctr.Final..rcv.glmnet Y: min < min of Train range: 1"
##      UniqueID Pplr.fctr.Final..rcv.glmnet WordCount.log1p WordCount.root2
## 8217     8217                           Y        1.609438               2
##                              id     cor.y exclude.as.feat cor.y.abs
## WordCount.log1p WordCount.log1p 0.2543196           FALSE 0.2543196
## WordCount.root2 WordCount.root2 0.2921207           FALSE 0.2921207
##                      cor.high.X freqRatio percentUnique zeroVar   nzv
## WordCount.log1p WordCount.root2  2.315789      24.15799   FALSE FALSE
## WordCount.root2            <NA>  2.315789      24.15799   FALSE FALSE
##                 is.cor.y.abs.low interaction.feat shapiro.test.p.value
## WordCount.log1p            FALSE               NA         1.576866e-49
## WordCount.root2            FALSE               NA         4.556481e-30
##                 rsp_var_raw id_var rsp_var       max min max.Pplr.fctr.N
## WordCount.log1p       FALSE     NA      NA   9.29771   0        8.819665
## WordCount.root2       FALSE     NA      NA 104.46052   0       82.249620
##                 max.Pplr.fctr.Y min.Pplr.fctr.N min.Pplr.fctr.Y
## WordCount.log1p         9.29771               0         1.94591
## WordCount.root2       104.46052               0         2.44949
##                 max.Pplr.fctr.All.X..rcv.glmnet.N
## WordCount.log1p                          7.069023
## WordCount.root2                         34.263683
##                 max.Pplr.fctr.All.X..rcv.glmnet.Y
## WordCount.log1p                          9.140883
## WordCount.root2                         96.581572
##                 min.Pplr.fctr.All.X..rcv.glmnet.N
## WordCount.log1p                                 0
## WordCount.root2                                 0
##                 min.Pplr.fctr.All.X..rcv.glmnet.Y
## WordCount.log1p                                 0
## WordCount.root2                                 0
##                 max.Pplr.fctr.Final..rcv.glmnet.N
## WordCount.log1p                           7.94094
## WordCount.root2                          53.00000
##                 max.Pplr.fctr.Final..rcv.glmnet.Y
## WordCount.log1p                          8.692322
## WordCount.root2                         77.175126
##                 min.Pplr.fctr.Final..rcv.glmnet.N
## WordCount.log1p                                 0
## WordCount.root2                                 0
##                 min.Pplr.fctr.Final..rcv.glmnet.Y
## WordCount.log1p                          1.609438
## WordCount.root2                          2.000000
## [1] "newobs Pplr.fctr.Final..rcv.glmnet Y: max > max of Train range: 684"
##      UniqueID Pplr.fctr.Final..rcv.glmnet PubDate.day.minutes.poly.1
## 6534     6534                           Y                 0.02047428
## 6535     6535                           Y                 0.02043797
## 6536     6536                           Y                 0.01840446
## 6537     6537                           Y                 0.01437377
## 6538     6538                           Y                 0.01408327
## 6539     6539                           Y                 0.01390170
##      PubDate.day.minutes.poly.3 PubDate.day.minutes.poly.4
## 6534                0.020825101                0.010136126
## 6535                0.020614888                0.009828175
## 6536                0.010192448               -0.003724687
## 6537               -0.003371562               -0.013710045
## 6538               -0.004024505               -0.013788116
## 6539               -0.004412350               -0.013803062
##      PubDate.day.minutes.poly.5 PubDate.juliandate PubDate.last32.log1p
## 6534               -0.004924177                335            10.036094
## 6535               -0.005269034                335            10.053458
## 6536               -0.017201046                335             9.939434
## 6537               -0.013755146                335             9.564863
## 6538               -0.012806242                335             9.542733
## 6539               -0.012191759                335             9.572898
##      WordCount.nexp
## 6534  4.609768e-243
## 6535   0.000000e+00
## 6536   0.000000e+00
## 6537   3.128062e-93
## 6538   0.000000e+00
## 6539   0.000000e+00
##      UniqueID Pplr.fctr.Final..rcv.glmnet PubDate.day.minutes.poly.1
## 6671     6671                           Y                0.002463239
## 6873     6873                           Y               -0.004799277
## 7317     7317                           Y                0.012231324
## 7528     7528                           Y                0.020437967
## 7734     7734                           Y                0.007002312
## 8145     8145                           Y               -0.010355102
##      PubDate.day.minutes.poly.3 PubDate.day.minutes.poly.4
## 6671               -0.005335173                0.006293548
## 6873                0.005333956                0.007895719
## 7317               -0.007282868               -0.012874650
## 7528                0.020614888                0.009828175
## 7734               -0.009404459               -0.002840648
## 8145                0.010795979               -0.003014347
##      PubDate.day.minutes.poly.5 PubDate.juliandate PubDate.last32.log1p
## 6671                0.007723689                336             9.941120
## 6873               -0.008268638                338            10.918446
## 7317               -0.006115494                345             9.800014
## 7528               -0.005269034                349            10.284797
## 7734                0.008456902                351             9.242904
## 8145               -0.009591936                357            11.002700
##      WordCount.nexp
## 6671   0.000000e+00
## 6873  2.032231e-313
## 7317  2.329036e-211
## 7528  2.699143e-152
## 7734  5.709040e-171
## 8145   0.000000e+00
##      UniqueID Pplr.fctr.Final..rcv.glmnet PubDate.day.minutes.poly.1
## 8386     8386                           Y               -0.004762964
## 8391     8391                           Y               -0.007559033
## 8392     8392                           Y               -0.007885846
## 8397     8397                           Y               -0.012134418
## 8399     8399                           Y               -0.014204235
## 8402     8402                           Y               -0.027458327
##      PubDate.day.minutes.poly.3 PubDate.day.minutes.poly.4
## 8386                0.005283252                0.007937277
## 8391                0.008724731                0.003454984
## 8392                0.009050092                0.002780751
## 8397                0.011101808               -0.007579609
## 8399                0.010203360               -0.012625915
## 8402               -0.044820236                0.033658267
##      PubDate.day.minutes.poly.5 PubDate.juliandate PubDate.last32.log1p
## 8386               -0.008201720                365             11.19043
## 8391               -0.011437989                365             11.15722
## 8392               -0.011513182                365             11.15937
## 8397               -0.005346258                365             11.11095
## 8399                0.002370447                365             11.10274
## 8402               -0.023627798                365             10.79296
##      WordCount.nexp
## 8386              0
## 8391              0
## 8392              0
## 8397              0
## 8399              0
## 8402              0
##                                                    id        cor.y
## PubDate.day.minutes.poly.1 PubDate.day.minutes.poly.1  0.156753478
## PubDate.day.minutes.poly.3 PubDate.day.minutes.poly.3  0.027983551
## PubDate.day.minutes.poly.4 PubDate.day.minutes.poly.4  0.073941394
## PubDate.day.minutes.poly.5 PubDate.day.minutes.poly.5 -0.055929231
## PubDate.juliandate                 PubDate.juliandate  0.014361075
## PubDate.last32.log1p             PubDate.last32.log1p  0.003558081
## WordCount.nexp                         WordCount.nexp -0.053208396
##                            exclude.as.feat   cor.y.abs         cor.high.X
## PubDate.day.minutes.poly.1           FALSE 0.156753478               <NA>
## PubDate.day.minutes.poly.3           FALSE 0.027983551               <NA>
## PubDate.day.minutes.poly.4           FALSE 0.073941394               <NA>
## PubDate.day.minutes.poly.5           FALSE 0.055929231               <NA>
## PubDate.juliandate                   FALSE 0.014361075 PubDate.month.fctr
## PubDate.last32.log1p                 FALSE 0.003558081               <NA>
## WordCount.nexp                       FALSE 0.053208396               <NA>
##                            freqRatio percentUnique zeroVar   nzv
## PubDate.day.minutes.poly.1   1.22549     18.080220   FALSE FALSE
## PubDate.day.minutes.poly.3   1.22549     18.080220   FALSE FALSE
## PubDate.day.minutes.poly.4   1.22549     18.080220   FALSE FALSE
## PubDate.day.minutes.poly.5   1.22549     18.080220   FALSE FALSE
## PubDate.juliandate           1.03252      1.393141   FALSE FALSE
## PubDate.last32.log1p         8.00000     90.998163   FALSE FALSE
## WordCount.nexp              17.76136     11.328843   FALSE FALSE
##                            is.cor.y.abs.low interaction.feat
## PubDate.day.minutes.poly.1            FALSE               NA
## PubDate.day.minutes.poly.3            FALSE               NA
## PubDate.day.minutes.poly.4            FALSE               NA
## PubDate.day.minutes.poly.5            FALSE               NA
## PubDate.juliandate                    FALSE               NA
## PubDate.last32.log1p                   TRUE               NA
## WordCount.nexp                        FALSE               NA
##                            shapiro.test.p.value rsp_var_raw id_var rsp_var
## PubDate.day.minutes.poly.1         1.590362e-18       FALSE     NA      NA
## PubDate.day.minutes.poly.3         9.822405e-52       FALSE     NA      NA
## PubDate.day.minutes.poly.4         1.523136e-47       FALSE     NA      NA
## PubDate.day.minutes.poly.5         1.157500e-41       FALSE     NA      NA
## PubDate.juliandate                 1.389406e-35       FALSE     NA      NA
## PubDate.last32.log1p               2.783236e-77       FALSE     NA      NA
## WordCount.nexp                     9.108805e-94       FALSE     NA      NA
##                                     max          min max.Pplr.fctr.N
## PubDate.day.minutes.poly.1   0.02475916  -0.02749464      0.02472285
## PubDate.day.minutes.poly.3   0.05215301  -0.04512497      0.05182988
## PubDate.day.minutes.poly.4   0.06677441  -0.01832740      0.06610094
## PubDate.day.minutes.poly.5   0.08471756  -0.02450918      0.08344228
## PubDate.juliandate         365.00000000 244.00000000    334.00000000
## PubDate.last32.log1p        12.32340669   0.00000000     12.21244232
## WordCount.nexp               1.00000000   0.00000000      1.00000000
##                            max.Pplr.fctr.Y min.Pplr.fctr.N min.Pplr.fctr.Y
## PubDate.day.minutes.poly.1    2.446866e-02     -0.02749464     -0.02745833
## PubDate.day.minutes.poly.3    4.959703e-02     -0.04512497     -0.04482024
## PubDate.day.minutes.poly.4    6.149053e-02     -0.01832740     -0.01821959
## PubDate.day.minutes.poly.5    7.481472e-02     -0.02450918     -0.02362780
## PubDate.juliandate            3.340000e+02    244.00000000    244.00000000
## PubDate.last32.log1p          1.221791e+01      0.00000000      0.00000000
## WordCount.nexp                2.478752e-03      0.00000000      0.00000000
##                            max.Pplr.fctr.All.X..rcv.glmnet.N
## PubDate.day.minutes.poly.1                        0.02450498
## PubDate.day.minutes.poly.3                        0.04991290
## PubDate.day.minutes.poly.4                        0.06213811
## PubDate.day.minutes.poly.5                        0.07601554
## PubDate.juliandate                              334.00000000
## PubDate.last32.log1p                             12.21791228
## WordCount.nexp                                    1.00000000
##                            max.Pplr.fctr.All.X..rcv.glmnet.Y
## PubDate.day.minutes.poly.1                        0.02472285
## PubDate.day.minutes.poly.3                        0.05182988
## PubDate.day.minutes.poly.4                        0.06610094
## PubDate.day.minutes.poly.5                        0.08344228
## PubDate.juliandate                              334.00000000
## PubDate.last32.log1p                             12.19622431
## WordCount.nexp                                    1.00000000
##                            min.Pplr.fctr.All.X..rcv.glmnet.N
## PubDate.day.minutes.poly.1                       -0.02749464
## PubDate.day.minutes.poly.3                       -0.04512497
## PubDate.day.minutes.poly.4                       -0.01832685
## PubDate.day.minutes.poly.5                       -0.02450918
## PubDate.juliandate                              244.00000000
## PubDate.last32.log1p                              0.00000000
## WordCount.nexp                                    0.00000000
##                            min.Pplr.fctr.All.X..rcv.glmnet.Y
## PubDate.day.minutes.poly.1                       -0.02745833
## PubDate.day.minutes.poly.3                       -0.04482024
## PubDate.day.minutes.poly.4                       -0.01791547
## PubDate.day.minutes.poly.5                       -0.02362780
## PubDate.juliandate                              244.00000000
## PubDate.last32.log1p                              0.00000000
## WordCount.nexp                                    0.00000000
##                            max.Pplr.fctr.Final..rcv.glmnet.N
## PubDate.day.minutes.poly.1                        0.02432341
## PubDate.day.minutes.poly.3                        0.04834381
## PubDate.day.minutes.poly.4                        0.05893666
## PubDate.day.minutes.poly.5                        0.07011504
## PubDate.juliandate                              365.00000000
## PubDate.last32.log1p                             12.32340669
## WordCount.nexp                                    1.00000000
##                            max.Pplr.fctr.Final..rcv.glmnet.Y
## PubDate.day.minutes.poly.1                        0.02475916
## PubDate.day.minutes.poly.3                        0.05215301
## PubDate.day.minutes.poly.4                        0.06677441
## PubDate.day.minutes.poly.5                        0.08471756
## PubDate.juliandate                              365.00000000
## PubDate.last32.log1p                             12.30546086
## WordCount.nexp                                    0.01831564
##                            min.Pplr.fctr.Final..rcv.glmnet.N
## PubDate.day.minutes.poly.1                       -0.02749464
## PubDate.day.minutes.poly.3                       -0.04512497
## PubDate.day.minutes.poly.4                       -0.01832268
## PubDate.day.minutes.poly.5                       -0.02450918
## PubDate.juliandate                              335.00000000
## PubDate.last32.log1p                              8.86290830
## WordCount.nexp                                    0.00000000
##                            min.Pplr.fctr.Final..rcv.glmnet.Y
## PubDate.day.minutes.poly.1                       -0.02745833
## PubDate.day.minutes.poly.3                       -0.04482024
## PubDate.day.minutes.poly.4                       -0.01820339
## PubDate.day.minutes.poly.5                       -0.02362780
## PubDate.juliandate                              335.00000000
## PubDate.last32.log1p                              8.83579237
## WordCount.nexp                                    0.00000000
## [1] "newobs total range outliers: 1870"
```

```
## numeric(0)
```

```
## [1] "glb_sel_mdl_id: All.X##rcv#glmnet"
```

```
## [1] "glb_fin_mdl_id: Final##rcv#glmnet"
```

```
## [1] "Cross Validation issues:"
##        MFO###myMFO_classfr  Random###myrandom_classfr 
##                          0                          0 
## Max.cor.Y.rcv.1X1###glmnet 
##                          0
```

```
##                                 max.Accuracy.OOB max.AUCROCR.OOB
## Max.cor.Y##rcv#rpart                   0.8200231       0.5892132
## Max.cor.Y.Time.Poly##rcv#glmnet        0.7754630       0.8049472
## All.X##rcv#glmnet                      0.7743056       0.7995951
## Interact.High.cor.Y##rcv#glmnet        0.7656250       0.8140971
## Max.cor.Y.rcv.1X1###glmnet             0.7604167       0.8116126
## Low.cor.X##rcv#glmnet                  0.7557870       0.8160359
## Max.cor.Y.Time.Lag##rcv#glmnet         0.6464120       0.8117635
## MFO###myMFO_classfr                    0.1331019       0.5000000
## Random###myrandom_classfr              0.1331019       0.4857956
## Final##rcv#glmnet                             NA              NA
##                                 max.AUCpROC.OOB max.Accuracy.fit
## Max.cor.Y##rcv#rpart                  0.5870523        0.9296422
## Max.cor.Y.Time.Poly##rcv#glmnet       0.5965780        0.9322790
## All.X##rcv#glmnet                     0.5950717        0.9320707
## Interact.High.cor.Y##rcv#glmnet       0.6009259        0.9315850
## Max.cor.Y.rcv.1X1###glmnet            0.5962443        0.9329725
## Low.cor.X##rcv#glmnet                 0.6015934        0.9323486
## Max.cor.Y.Time.Lag##rcv#glmnet        0.5925379        0.9279084
## MFO###myMFO_classfr                   0.5000000        0.1796420
## Random###myrandom_classfr             0.5125675        0.1796420
## Final##rcv#glmnet                            NA        0.9065624
##                                 opt.prob.threshold.fit
## Max.cor.Y##rcv#rpart                               0.6
## Max.cor.Y.Time.Poly##rcv#glmnet                    0.4
## All.X##rcv#glmnet                                  0.4
## Interact.High.cor.Y##rcv#glmnet                    0.4
## Max.cor.Y.rcv.1X1###glmnet                         0.5
## Low.cor.X##rcv#glmnet                              0.2
## Max.cor.Y.Time.Lag##rcv#glmnet                     0.3
## MFO###myMFO_classfr                                0.1
## Random###myrandom_classfr                          0.1
## Final##rcv#glmnet                                  0.2
##                                 opt.prob.threshold.OOB
## Max.cor.Y##rcv#rpart                               0.6
## Max.cor.Y.Time.Poly##rcv#glmnet                    0.1
## All.X##rcv#glmnet                                  0.1
## Interact.High.cor.Y##rcv#glmnet                    0.1
## Max.cor.Y.rcv.1X1###glmnet                         0.1
## Low.cor.X##rcv#glmnet                              0.1
## Max.cor.Y.Time.Lag##rcv#glmnet                     0.1
## MFO###myMFO_classfr                                0.1
## Random###myrandom_classfr                          0.1
## Final##rcv#glmnet                                   NA
```

```
## [1] "All.X##rcv#glmnet OOB confusion matrix & accuracy: "
##          Prediction
## Reference    N    Y
##         N 1195  303
##         Y   87  143
##                         .freqRatio.Fit .freqRatio.OOB .freqRatio.Tst
## OpEd#Opnn#                 0.090965862   0.0515046296    0.087700535
## Styls#U.S.#                0.026436303   0.0289351852    0.032620321
## Scnc#Hlth#                 0.030807660   0.0277777778    0.030481283
## #Opnn#ThPblcEdtr           0.003330558   0.0023148148    0.005347594
## Bsnss#Crsswrds/Gms#        0.021856786   0.0104166667    0.022459893
## Bsnss#Tchnlgy#             0.044338052   0.0729166667    0.060962567
## #Opnn#RmFrDbt              0.008742714   0.0115740741    0.010695187
## Bsnss#BsnssDy#Dlbk         0.130932556   0.1869212963    0.162566845
## ##                         0.190049958   0.2146990741    0.182887701
## Cltr#Arts#                 0.101998335   0.1070601852    0.093048128
## Mtr#N.Y./Rgn#              0.026644463   0.0405092593    0.035828877
## Styls##Fshn                0.021648626   0.0086805556    0.008021390
## Bsnss#BsnssDy#SmllBsnss    0.020815987   0.0231481481    0.021925134
## Frgn#Wrld#AsPcfc           0.031223980   0.0306712963    0.029946524
## TStyl##                    0.129683597   0.0584490741    0.056149733
## Trvl#Trvl#                 0.017277269   0.0196759259    0.018716578
## #Mltmd#                    0.019150708   0.0283564815    0.027807487
## myOthr                     0.006869276   0.0028935185    0.002673797
## #U.S.#Edctn                0.050582848   0.0474537037    0.047593583
## Cltr##                              NA   0.0005787037    0.037433155
## Frgn#Wrld#                 0.026644463   0.0254629630    0.025133690
##                         .n.Fit .n.New.N .n.New.Y .n.OOB .n.Trn.N .n.Trn.Y
## OpEd#Opnn#                 437       NA      164     89      117      409
## Styls#U.S.#                127       NA       61     50       77      100
## Scnc#Hlth#                 148       NA       57     48       74      122
## #Opnn#ThPblcEdtr            16       NA       10      4        4       16
## Bsnss#Crsswrds/Gms#        105       NA       42     18       20      103
## Bsnss#Tchnlgy#             213       34       80    126      288       51
## #Opnn#RmFrDbt               42       19        1     20       61        1
## Bsnss#BsnssDy#Dlbk         629      197      107    323      864       88
## ##                         913      265       77    371     1169      115
## Cltr#Arts#                 490      157       17    185      625       50
## Mtr#N.Y./Rgn#              128       37       30     70      181       17
## Styls##Fshn                104       15       NA     15      118        1
## Bsnss#BsnssDy#SmllBsnss    100       35        6     40      135        5
## Frgn#Wrld#AsPcfc           150       49        7     53      200        3
## TStyl##                    623      104        1    101      715        9
## Trvl#Trvl#                  83       35       NA     34      116        1
## #Mltmd#                     92       52       NA     49      139        2
## myOthr                      33        5       NA      5       38       NA
## #U.S.#Edctn                243       87        2     82      325       NA
## Cltr##                      NA       48       22      1        1       NA
## Frgn#Wrld#                 128       47       NA     44      172       NA
##                         .n.Tst .n.fit .n.new .n.trn err.abs.OOB.mean
## OpEd#Opnn#                 164    437    164    526       0.57506350
## Styls#U.S.#                 61    127     61    177       0.51995217
## Scnc#Hlth#                  57    148     57    196       0.51438534
## #Opnn#ThPblcEdtr            10     16     10     20       0.49777421
## Bsnss#Crsswrds/Gms#         42    105     42    123       0.48547153
## Bsnss#Tchnlgy#             114    213    114    339       0.20447524
## #Opnn#RmFrDbt               20     42     20     62       0.19588822
## Bsnss#BsnssDy#Dlbk         304    629    304    952       0.18546652
## ##                         342    913    342   1284       0.18256154
## Cltr#Arts#                 174    490    174    675       0.16900564
## Mtr#N.Y./Rgn#               67    128     67    198       0.16276609
## Styls##Fshn                 15    104     15    119       0.12795213
## Bsnss#BsnssDy#SmllBsnss     41    100     41    140       0.11571966
## Frgn#Wrld#AsPcfc            56    150     56    203       0.09046102
## TStyl##                    105    623    105    724       0.07789865
## Trvl#Trvl#                  35     83     35    117       0.07507486
## #Mltmd#                     52     92     52    141       0.07179209
## myOthr                       5     33      5     38       0.06786732
## #U.S.#Edctn                 89    243     89    325       0.06237106
## Cltr##                      70     NA     70      1       0.04526133
## Frgn#Wrld#                  47    128     47    172       0.04276137
##                         err.abs.fit.mean err.abs.new.mean err.abs.trn.mean
## OpEd#Opnn#                    0.27422533               NA       0.33251375
## Styls#U.S.#                   0.44589721               NA       0.45567870
## Scnc#Hlth#                    0.38524093               NA       0.38793668
## #Opnn#ThPblcEdtr              0.41155588               NA       0.32998288
## Bsnss#Crsswrds/Gms#           0.25730899               NA       0.23852271
## Bsnss#Tchnlgy#                0.18719813               NA       0.23001139
## #Opnn#RmFrDbt                 0.16749556               NA       0.06272545
## Bsnss#BsnssDy#Dlbk            0.13007797               NA       0.15397456
## ##                            0.11064396               NA       0.13393000
## Cltr#Arts#                    0.09840114               NA       0.11214813
## Mtr#N.Y./Rgn#                 0.12031512               NA       0.14888110
## Styls##Fshn                   0.06526513               NA       0.03797197
## Bsnss#BsnssDy#SmllBsnss       0.10499902               NA       0.08673011
## Frgn#Wrld#AsPcfc              0.09269842               NA       0.04536221
## TStyl##                       0.05361215               NA       0.02965505
## Trvl#Trvl#                    0.04822508               NA       0.03739970
## #Mltmd#                       0.06404709               NA       0.04594878
## myOthr                        0.07031647               NA       0.03739458
## #U.S.#Edctn                   0.05014178               NA       0.01357338
## Cltr##                                NA               NA       0.06664375
## Frgn#Wrld#                    0.04253820               NA       0.01777674
##                         err.abs.OOB.sum err.abs.fit.sum err.abs.new.sum
## OpEd#Opnn#                  51.18065179      119.836470              NA
## Styls#U.S.#                 25.99760848       56.628946              NA
## Scnc#Hlth#                  24.69049636       57.015658              NA
## #Opnn#ThPblcEdtr             1.99109686        6.584894              NA
## Bsnss#Crsswrds/Gms#          8.73848753       27.017444              NA
## Bsnss#Tchnlgy#              25.76388022       39.873202              NA
## #Opnn#RmFrDbt                3.91776433        7.034814              NA
## Bsnss#BsnssDy#Dlbk          59.90568728       81.819044              NA
## ##                          67.73033042      101.017935              NA
## Cltr#Arts#                  31.26604266       48.216560              NA
## Mtr#N.Y./Rgn#               11.39362621       15.400336              NA
## Styls##Fshn                  1.91928194        6.787573              NA
## Bsnss#BsnssDy#SmllBsnss      4.62878656       10.499902              NA
## Frgn#Wrld#AsPcfc             4.79443402       13.904763              NA
## TStyl##                      7.86776348       33.400369              NA
## Trvl#Trvl#                   2.55254529        4.002682              NA
## #Mltmd#                      3.51781240        5.892332              NA
## myOthr                       0.33933658        2.320443              NA
## #U.S.#Edctn                  5.11442705       12.184453              NA
## Cltr##                       0.04526133              NA              NA
## Frgn#Wrld#                   1.88150036        5.444890              NA
##                         err.abs.trn.sum
## OpEd#Opnn#                 174.90223507
## Styls#U.S.#                 80.65513073
## Scnc#Hlth#                  76.03558884
## #Opnn#ThPblcEdtr             6.59965754
## Bsnss#Crsswrds/Gms#         29.33829341
## Bsnss#Tchnlgy#              77.97386000
## #Opnn#RmFrDbt                3.88897815
## Bsnss#BsnssDy#Dlbk         146.58378254
## ##                         171.96611762
## Cltr#Arts#                  75.69998885
## Mtr#N.Y./Rgn#               29.47845712
## Styls##Fshn                  4.51866405
## Bsnss#BsnssDy#SmllBsnss     12.14221607
## Frgn#Wrld#AsPcfc             9.20852860
## TStyl##                     21.47025749
## Trvl#Trvl#                   4.37576471
## #Mltmd#                      6.47877859
## myOthr                       1.42099408
## #U.S.#Edctn                  4.41134846
## Cltr##                       0.06664375
## Frgn#Wrld#                   3.05759995
##   .freqRatio.Fit   .freqRatio.OOB   .freqRatio.Tst           .n.Fit 
##               NA         1.000000         1.000000               NA 
##         .n.New.N         .n.New.Y           .n.OOB         .n.Trn.N 
##               NA               NA      1728.000000      5439.000000 
##         .n.Trn.Y           .n.Tst           .n.fit           .n.new 
##               NA      1870.000000               NA      1870.000000 
##           .n.trn err.abs.OOB.mean err.abs.fit.mean err.abs.new.mean 
##      6532.000000         4.469969               NA               NA 
## err.abs.trn.mean  err.abs.OOB.sum  err.abs.fit.sum  err.abs.new.sum 
##         3.004762       345.236821               NA               NA 
##  err.abs.trn.sum 
##       940.272886
```

```
##                                     All.X__rcv_glmnet.imp
## PubDate.day.minutes.poly.1                     100.000000
## NDSS.my.fctrOpEd#Opnn#                          68.562523
## NDSS.my.fctrBsnss#Crsswrds/Gms#                 63.559761
## NDSS.my.fctrScnc#Hlth#                          54.307908
## NDSS.my.fctrStyls#U.S.#                         50.421526
## NDSS.my.fctr#Opnn#ThPblcEdtr                    48.417943
## WordCount.log1p                                  6.701272
## WordCount.root2                                  6.021636
## NDSS.my.fctrBsnss#Tchnlgy#                       5.690755
## PubDate.wkend                                    5.345889
## PubDate.day.minutes.poly.2                       5.307609
## PubDate.day.minutes.poly.4                       5.307609
## PubDate.day.minutes.poly.3                       5.307609
## PubDate.hour.fctr(7.67,15.3]                     5.307609
## PubDate.last16.log1p                             5.307609
## PubDate.wkday.fctr1                              5.307609
## PubDate.last32.log1p                             5.307609
## PubDate.date.fctr(7,13]                          5.307609
## .rnorm                                           5.307609
## NDSS.my.fctrCltr##                               5.307609
## NDSS.my.fctrMtr#N.Y./Rgn#                        5.307609
## PubDate.date.fctr(19,25]                         5.307609
## PubDate.date.fctr(25,31]                         5.307609
## PubDate.day.minutes.poly.5                       5.307609
## PubDate.hour.fctr(15.3,23]                       5.307609
## PubDate.juliandate                               5.307609
## PubDate.last2.log1p                              5.307609
## PubDate.last4.log1p                              5.307609
## PubDate.last8.log1p                              5.307609
## PubDate.minute.fctr(14.8,29.5]                   5.307609
## PubDate.minute.fctr(44.2,59.1]                   5.307609
## PubDate.month.fctr10                             5.307609
## PubDate.month.fctr12                             5.307609
## PubDate.second.fctr(14.8,29.5]                   5.307609
## PubDate.second.fctr(29.5,44.2]                   5.307609
## PubDate.wkday.fctr3                              5.307609
## PubDate.wkday.fctr4                              5.307609
## WordCount.nexp                                   5.307609
## PubDate.wkday.fctr2                              5.307609
## PubDate.date.fctr(13,19]                         5.307609
## NDSS.my.fctrCltr#Arts#                           5.307609
## PubDate.month.fctr11                             5.307609
## PubDate.second.fctr(44.2,59.1]                   5.307609
## PubDate.wkday.fctr6                              5.307609
## PubDate.wkday.fctr5                              5.307609
## PubDate.minute.fctr(29.5,44.2]                   5.307609
## NDSS.my.fctrBsnss#BsnssDy#Dlbk                   5.307609
## NDSS.my.fctrTrvl#Trvl#                           5.307609
## NDSS.my.fctrBsnss#BsnssDy#SmllBsnss              5.307609
## NDSS.my.fctr#Mltmd#                              5.307609
##                                     Final__rcv_glmnet.imp
## PubDate.day.minutes.poly.1                       90.40903
## NDSS.my.fctrOpEd#Opnn#                           33.93385
## NDSS.my.fctrBsnss#Crsswrds/Gms#                  31.44233
## NDSS.my.fctrScnc#Hlth#                           28.67549
## NDSS.my.fctrStyls#U.S.#                          27.49488
## NDSS.my.fctr#Opnn#ThPblcEdtr                     30.81095
## WordCount.log1p                                  16.43344
## WordCount.root2                                  14.82783
## NDSS.my.fctrBsnss#Tchnlgy#                       17.06756
## PubDate.wkend                                    16.18021
## PubDate.day.minutes.poly.2                      100.00000
## PubDate.day.minutes.poly.4                       32.18318
## PubDate.day.minutes.poly.3                       25.72759
## PubDate.hour.fctr(7.67,15.3]                     15.42835
## PubDate.last16.log1p                             15.05236
## PubDate.wkday.fctr1                              15.02383
## PubDate.last32.log1p                             14.63405
## PubDate.date.fctr(7,13]                          14.56304
## .rnorm                                           14.55456
## NDSS.my.fctrCltr##                               14.55456
## NDSS.my.fctrMtr#N.Y./Rgn#                        14.55456
## PubDate.date.fctr(19,25]                         14.55456
## PubDate.date.fctr(25,31]                         14.55456
## PubDate.day.minutes.poly.5                       14.55456
## PubDate.hour.fctr(15.3,23]                       14.55456
## PubDate.juliandate                               14.55456
## PubDate.last2.log1p                              14.55456
## PubDate.last4.log1p                              14.55456
## PubDate.last8.log1p                              14.55456
## PubDate.minute.fctr(14.8,29.5]                   14.55456
## PubDate.minute.fctr(44.2,59.1]                   14.55456
## PubDate.month.fctr10                             14.55456
## PubDate.month.fctr12                             14.55456
## PubDate.second.fctr(14.8,29.5]                   14.55456
## PubDate.second.fctr(29.5,44.2]                   14.55456
## PubDate.wkday.fctr3                              14.55456
## PubDate.wkday.fctr4                              14.55456
## WordCount.nexp                                   14.55456
## PubDate.wkday.fctr2                              14.50995
## PubDate.date.fctr(13,19]                         14.31658
## NDSS.my.fctrCltr#Arts#                           14.24518
## PubDate.month.fctr11                             14.24288
## PubDate.second.fctr(44.2,59.1]                   14.23286
## PubDate.wkday.fctr6                              14.02941
## PubDate.wkday.fctr5                              13.90528
## PubDate.minute.fctr(29.5,44.2]                   13.84506
## NDSS.my.fctrBsnss#BsnssDy#Dlbk                   13.73612
## NDSS.my.fctrTrvl#Trvl#                           12.16166
## NDSS.my.fctrBsnss#BsnssDy#SmllBsnss              11.40620
## NDSS.my.fctr#Mltmd#                              10.67411
```

```
## [1] "glbObsNew prediction stats:"
```

```
## 
##    N    Y 
## 1186  684
```

```
##                   label step_major step_minor label_minor     bgn     end
## 16     predict.data.new          8          0           0 269.154 285.003
## 17 display.session.info          9          0           0 285.004      NA
##    elapsed
## 16   15.85
## 17      NA
```

Null Hypothesis ($\sf{H_{0}}$): mpg is not impacted by am_fctr.  
The variance by am_fctr appears to be independent. 
#```{r q1, cache=FALSE}
# print(t.test(subset(cars_df, am_fctr == "automatic")$mpg, 
#              subset(cars_df, am_fctr == "manual")$mpg, 
#              var.equal=FALSE)$conf)
#```
We reject the null hypothesis i.e. we have evidence to conclude that am_fctr impacts mpg (95% confidence). Manual transmission is better for miles per gallon versus automatic transmission.


```
##                      label step_major step_minor label_minor     bgn
## 10              fit.models          6          0           0  80.360
## 14       fit.data.training          7          0           0 192.948
## 11              fit.models          6          1           1 150.660
## 5         extract.features          3          0           0  26.626
## 9          select.features          5          0           0  57.659
## 16        predict.data.new          8          0           0 269.154
## 12              fit.models          6          2           2 174.229
## 15       fit.data.training          7          1           1 258.677
## 1              import.data          1          0           0  11.375
## 6      manage.missing.data          3          1           1  49.558
## 13              fit.models          6          3           3 187.078
## 2             inspect.data          2          0           0  21.011
## 8  partition.data.training          4          0           0  56.152
## 3               scrub.data          2          1           1  25.368
## 4           transform.data          2          2           2  26.529
## 7             cluster.data          3          2           2  56.087
##        end elapsed duration
## 10 150.659  70.299   70.299
## 14 258.677  65.729   65.729
## 11 174.228  23.568   23.568
## 5   49.557  22.932   22.931
## 9   80.360  22.701   22.701
## 16 285.003  15.850   15.849
## 12 187.078  12.849   12.849
## 15 269.153  10.476   10.476
## 1   21.011   9.636    9.636
## 6   56.087   6.529    6.529
## 13 192.947   5.869    5.869
## 2   25.368   4.357    4.357
## 8   57.658   1.507    1.506
## 3   26.528   1.160    1.160
## 4   26.626   0.097    0.097
## 7   56.151   0.064    0.064
## [1] "Total Elapsed Time: 285.003 secs"
```

![](NYTBlogs3_feat_PubDate_files/figure-html/display.session.info-1.png) 

```
##                              label step_major step_minor      label_minor
## 8           fit.models_0_Low.cor.X          1          7           glmnet
## 4   fit.models_0_Max.cor.Y.rcv.*X*          1          3           glmnet
## 7 fit.models_0_Interact.High.cor.Y          1          6           glmnet
## 5 fit.models_0_Max.cor.Y.Time.Poly          1          4           glmnet
## 6  fit.models_0_Max.cor.Y.Time.Lag          1          5           glmnet
## 3              fit.models_0_Random          1          2 myrandom_classfr
## 2                 fit.models_0_MFO          1          1    myMFO_classfr
## 1                 fit.models_0_bgn          1          0            setup
##       bgn     end elapsed duration
## 8 136.861 150.646  13.785   13.785
## 4  88.902 102.249  13.347   13.347
## 7 124.990 136.861  11.871   11.871
## 5 102.250 114.112  11.863   11.862
## 6 114.113 124.990  10.877   10.877
## 3  84.532  88.901   4.370    4.369
## 2  81.489  84.531   3.042    3.042
## 1  81.459  81.488   0.030    0.029
## [1] "Total Elapsed Time: 150.646 secs"
```

![](NYTBlogs3_feat_PubDate_files/figure-html/display.session.info-2.png) 
