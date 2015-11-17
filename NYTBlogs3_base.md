# NYTimes:Blogs: Popular classification:: NYTBlogs3_base
bdanalytics  

**  **    
**Date: (Mon) Nov 16, 2015**    

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
- inspect.data chunk:
    - For date variables
        - Appropriate factors ?
        - Different / More last* features ?
        
- scrub.data chunk:        
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
    - move model_type parameter to myfit_mdl before indep_vars_vctr (keep all model_* together)
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
- Replicate myfit_mdl_classification features in myfit_mdl_regression
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
glb_rsp_var <- "Popular.fctr" # glb_rsp_var_raw # or "Popular.fctr"

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
glb_map_rsp_raw_to_var(tst <- c(NA, 0, 1)) 
```

```
## [1] <NA> N    Y   
## Levels: N Y
```

```r
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
glb_map_rsp_var_to_raw(glb_map_rsp_raw_to_var(tst))
```

```
## [1] NA  0  1
```

```r
if ((glb_rsp_var != glb_rsp_var_raw) && is.null(glb_map_rsp_raw_to_var))
    stop("glb_map_rsp_raw_to_var function expected")
glb_rsp_var_out <- paste0(glb_rsp_var, ".predict.") # mdl_id is appended later

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
    , "WordCount" # Feature Engineering done with prior features
    , "Headline", "Snippet", "Abstract", "PubDate"
                    ) 
if (glb_rsp_var_raw != glb_rsp_var)
    glbFeatsExclude <- union(glbFeatsExclude, glb_rsp_var_raw)                    

glbFeatsInteractionOnly <- list()
#glbFeatsInteractionOnly[["carrier.fctr"]] <- "cellular.fctr"

# currently does not handle more than 1 column; consider concatenating multiple columns
glb_id_var <- c(NULL 
                , "UniqueID")
glb_category_var <- "NDSSName.my.fctr" # choose from c(NULL : default, "<category>")

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

# Add logs of numerics that are not distributed normally
#   Derive & keep multiple transformations of the same feature, if normality is hard to achieve with just one transformation
#   Right skew: logp1; sqrt; ^ 1/3; logp1(logp1); log10; exp(-<feat>/constant)

# glbFeatsDerive[["<feat.my.sfx>"]] <- list(
    # character
#     mapfn = function(Week) { return(substr(Week, 1, 10)) }

#print(mod_raw <- grep("&#034;", glbObsAll[, txt_var], value = TRUE)) 
#print(mod_raw <- glbObsAll[c(88,187,280,1040,1098), txt_var])
#print(mod_raw <- glbObsAll[sel_obs(list(descr.my.contains="\\bdoes( +)not\\b")), glbFeatsText])
#print(mod_raw <- glbObsAll[sel_obs(list(descr.my.contains="\\bipad [[:digit:]]\\b")), glbFeatsText][01:10])
#print(mod_raw <- glbObsAll[sel_obs(list(descr.my.contains="pad mini")), glbFeatsText][11:20])
#print(mod_raw <- glbObsAll[sel_obs(list(descr.my.contains="pad mini")), glbFeatsText][21:30])
#print(mod_raw <- glbObsAll[sel_obs(list(descr.my.contains="pad mini")), glbFeatsText][31:40])
#glbObsAll[which(glb_post_stop_words_terms_mtrx_lst[[txt_var]][, subset(glb_post_stop_words_terms_df_lst[[txt_var]], term %in% c("conditionminimal"))$pos] > 0), "description"]

#     mapfn = function(description) { mod_raw <- description;

    # This is here because it does not work if it's in txt_map_filename
#         mod_raw <- gsub(paste0(c("\n", "\211", "\235", "\317", "\333"), collapse = "|"), " ", mod_raw)

    # Don't parse for "." because of ".com"; use customized gsub for that text
#         mod_raw <- gsub("(\\w)(!|\\*|,|-|/)(\\w)", "\\1\\2 \\3", mod_raw);

#         return(mod_raw) }

    # numeric
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
#     mapfn = function(descriptor) { return(plyr::revalue(descriptor, c(
#         "ABANDONED BUILDING"  = "OTHER",
#         "##"                  = "##"
#                                           ))) }
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
#         paste(gsub(" ", "", productline), as.numeric(nchar(description) > 0), sep = "#")) }

#     , args = c("<arg1>", "<arg2>"))

# # If glbObsAll is not sorted in the desired manner
#     mapfn=function(Week) { return(coredata(lag(zoo(orderBy(~Week, glbObsAll)$ILI), -2, na.pad=TRUE))) }
#     mapfn=function(ILI) { return(coredata(lag(zoo(ILI), -2, na.pad=TRUE))) }
#     mapfn=function(ILI.2.lag) { return(log(ILI.2.lag)) }

# glbFeatsDerive[["dummy.my"]] <- list(
#     mapfn = function(UniqueID) { return(UniqueID) }       
#     , args = c("UniqueID"))    
glbFeatsDerive[["NDSSName.my"]] <- list(
    mapfn = function(NewsDesk, SectionName, SubsectionName) { 
        descriptor <- 
            gsub(" ", "", paste(NewsDesk, SectionName, SubsectionName, sep = "#"))
        return(plyr::revalue(descriptor, c(NULL
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
                                          )))
        }
    , args = c("NewsDesk", "SectionName", "SubsectionName"))    

# glbFeatsDerive[["<feat.my.sfx>"]] <- list(
#     mapfn = function(<arg1>, <arg2>) { return(function(<arg1>, <arg2>)) } 
#   , args = c("<arg1>", "<arg2>"))

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

# glbFeatsDerive[["<var1>"]] <- glbFeatsDerive[["<var2>"]]

glb_derive_vars <- names(glbFeatsDerive)
# tst <- "descr.my"; args_lst <- NULL; for (arg in glbFeatsDerive[[tst]]$args) args_lst[[arg]] <- glbObsAll[, arg]; print(head(args_lst[[arg]])); print(head(drv_vals <- do.call(glbFeatsDerive[[tst]]$mapfn, args_lst))); 
# print(which_ix <- which(args_lst[[arg]] == 0.75)); print(drv_vals[which_ix]); 

glb_date_vars <- NULL # or c("<date_var>")
glb_date_fmts <- list(); #glb_date_fmts[["<date_var>"]] <- "%m/%e/%y"
glb_date_tzs <- list();  #glb_date_tzs[["<date_var>"]] <- "America/New_York"
#grep("America/New", OlsonNames(), value=TRUE)

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
#glbObsAll[which(terms_stop_mtrx[, grep("16", dimnames(terms_stop_mtrx)$Terms)[1]] > 0), c(glb_category_var, "storage", txt_var)]

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
#glbObsAll[which(TfIdf_stem_mtrx[, 191] > 0), c(glb_id_var, glb_category_var, txt_var)]
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
glb_impute_na_data <- FALSE # or TRUE
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

#mdlId <- "CSM.X.glm"; vars <- myextract_actual_feats(row.names(orderBy(reformulate(c("-", paste0(mdlId, ".importance"))), myget_feats_importance(glb_models_lst[[mdlId]])))); 
#model_diags_df <- glb_get_predictions(model_diags_df, mdlId,glb_rsp_var_out)
#obs_ix <- row.names(model_diags_df) %in% names(outliers$rstudent)[1]
#obs_ix <- which(is.na(model_diags_df$.rstudent))
#obs_ix <- which(is.na(model_diags_df$.dffits))
#myplot_parcoord(obs_df=model_diags_df[, c(glb_id_var, glb_category_var, ".rstudent", ".dffits", ".hatvalues", glb_rsp_var, paste0(glb_rsp_var_out, mdlId), vars[1:min(20, length(vars))])], obs_ix=obs_ix, id_var=glb_id_var, category_var=glb_category_var)

#model_diags_df[row.names(model_diags_df) %in% names(outliers$rstudent)[c(1:2)], ]
#ctgry_diags_df <- model_diags_df[model_diags_df[, glb_category_var] %in% c("Unknown#0"), ]
#myplot_parcoord(obs_df=ctgry_diags_df[, c(glb_id_var, glb_category_var, ".rstudent", ".dffits", ".hatvalues", glb_rsp_var, "startprice.log10.predict.RFE.X.glmnet", indep_vars[1:20])], obs_ix=row.names(ctgry_diags_df) %in% names(outliers$rstudent)[1], id_var=glb_id_var, category_var=glb_category_var)
#table(glbObsFit[model_diags_df[, glb_category_var] %in% c("iPad1#1"), "startprice.log10.cut.fctr"])
#glbObsFit[model_diags_df[, glb_category_var] %in% c("iPad1#1"), c(glb_id_var, "startprice")]

# No outliers & .dffits == NaN
#myplot_parcoord(obs_df=model_diags_df[, c(glb_id_var, glb_category_var, glb_rsp_var, "startprice.log10.predict.RFE.X.glmnet", indep_vars[1:10])], obs_ix=seq(1:nrow(model_diags_df))[is.na(model_diags_df$.dffits)], id_var=glb_id_var, category_var=glb_category_var)

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
glb_tune_models_df <- data.frame()

    #avNNet    
    #   size=[1] 3 5 7 9; decay=[0] 1e-04 0.001  0.01   0.1; bag=[FALSE]; RMSE=1.3300906 

    #bagEarth
    #   degree=1 [2] 3; nprune=64 128 256 512 [1024]; RMSE=0.6486663 (up)
# glb_tune_models_df <- myrbind_df(glb_tune_models_df, rbind(data.frame()
#     ,data.frame(method = "bagEarth", parameter = "nprune", vals = "256")
#     ,data.frame(method = "bagEarth", parameter = "degree", vals = "2")    
# ))

    #earth 
    #   degree=[1]; nprune=2  [9] 17 25 33; RMSE=0.1334478
    
    #gbm 
    #   shrinkage=0.05 [0.10] 0.15 0.20 0.25; n.trees=100 150 200 [250] 300; interaction.depth=[1] 2 3 4 5; n.minobsinnode=[10]; RMSE=0.2008313     
# glb_tune_models_df <- myrbind_df(glb_tune_models_df, rbind(data.frame()
#     ,data.frame(method = "gbm", parameter = "shrinkage", min = 0.05, max = 0.25, by = 0.05)
#     ,data.frame(method = "gbm", parameter = "n.trees", min = 100, max = 300, by = 50)
#     ,data.frame(method = "gbm", parameter = "interaction.depth", min = 1, max = 5, by = 1)
#     ,data.frame(method = "gbm", parameter = "n.minobsinnode", min = 10, max = 10, by = 10)
#     #seq(from=0.05,  to=0.25, by=0.05)
# ))

    #glmnet
    #   alpha=0.100 [0.325] 0.550 0.775 1.000; lambda=0.0005232693 0.0024288010 0.0112734954 [0.0523269304] 0.2428800957; RMSE=0.6164891
# glb_tune_models_df <- myrbind_df(glb_tune_models_df, rbind(data.frame()
#     ,data.frame(method = "glmnet", parameter = "alpha", vals = "0.550 0.775 0.8875 0.94375 1.000")
#     ,data.frame(method = "glmnet", parameter = "lambda", vals = "9.858855e-05 0.0001971771 0.0009152152 0.0042480525 0.0197177130")    
# ))

    #nnet    
    #   size=3 5 [7] 9 11; decay=0.0001 0.001 0.01 [0.1] 0.2; RMSE=0.9287422
# glb_tune_models_df <- myrbind_df(glb_tune_models_df, rbind(data.frame()
#     ,data.frame(method = "nnet", parameter = "size", vals = "3 5 7 9 11")
#     ,data.frame(method = "nnet", parameter = "decay", vals = "0.0001 0.0010 0.0100 0.1000 0.2000")    
# ))

    #rf # Don't bother; results are not deterministic
    #       mtry=2  35  68 [101] 134; RMSE=0.1339974
# glb_tune_models_df <- myrbind_df(glb_tune_models_df, rbind(data.frame()
#     ,data.frame(method = "rf", parameter = "mtry", vals = "2 5 9 13 17")
# ))

    #rpart 
    #   cp=0.020 [0.025] 0.030 0.035 0.040; RMSE=0.1770237
# glb_tune_models_df <- myrbind_df(glb_tune_models_df, rbind(data.frame()    
#     ,data.frame(method = "rpart", parameter = "cp", vals = "0.004347826 0.008695652 0.017391304 0.021739130 0.034782609")
# ))
    
    #svmLinear
    #   C=0.01 0.05 [0.10] 0.50 1.00 2.00 3.00 4.00; RMSE=0.1271318; 0.1296718
# glb_tune_models_df <- myrbind_df(glb_tune_models_df, rbind(data.frame()
#     ,data.frame(method = "svmLinear", parameter = "C", vals = "0.01 0.05 0.1 0.5 1")
# ))

    #svmLinear2    
    #   cost=0.0625 0.1250 [0.25] 0.50 1.00; RMSE=0.1276354 
# glb_tune_models_df <- myrbind_df(glb_tune_models_df, rbind(data.frame()
#     ,data.frame(method = "svmLinear2", parameter = "cost", vals = "0.0625 0.125 0.25 0.5 1")
# ))

    #svmPoly    
    #   degree=[1] 2 3 4 5; scale=0.01 0.05 [0.1] 0.5 1; C=0.50 1.00 [2.00] 3.00 4.00; RMSE=0.1276130
# glb_tune_models_df <- myrbind_df(glb_tune_models_df, rbind(data.frame()
#     ,data.frame(method="svmPoly", parameter="degree", min=1, max=5, by=1) #seq(1, 5, 1)
#     ,data.frame(method="svmPoly", parameter="scale", vals="0.01, 0.05, 0.1, 0.5, 1")
#     ,data.frame(method="svmPoly", parameter="C", vals="0.50, 1.00, 2.00, 3.00, 4.00")    
# ))

    #svmRadial
    #   sigma=[0.08674323]; C=0.25 0.50 1.00 [2.00] 4.00; RMSE=0.1614957
    
#glb_to_sav(); all.equal(sav_models_df, glb_models_df)
    
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
# tmp_fitobs_df <- glbObsFit[, grep(paste0("^", gsub(".", "\\.", glb_rsp_var_out, fixed = TRUE), "CSM\\.X\\.(.*)\\.prob"), names(glbObsFit), value = TRUE)]; cor_mtrx <- cor(tmp_fitobs_df); cor_vctr <- sort(cor_mtrx[row.names(orderBy(~-Overall, varImp(glb_models_lst[["Ensemble.repeatedcv.glmnet"]])$importance))[1], ]); summary(cor_vctr); cor_vctr
#ntv.glm <- glm(reformulate(indep_vars, glb_rsp_var), family = "binomial", data = glbObsFit)
#step.glm <- step(ntv.glm)

glb_sel_mdl_id <- "All.X.glmnet" #select from c(NULL, "All.X.glmnet", "RFE.X.glmnet", <mdlId>)
glb_fin_mdl_id <- NULL #select from c(NULL, glb_sel_mdl_id)

glb_dsp_cols <- c(NULL
#               List critical cols excl. glb_id_var, glb_category_var & glb_rsp_var
                  )

# Output specs
glbOutDataVizFname <- "NYTBlogs3_obsall.csv" # choose from c(NULL, "NYTBlogs3_obsall.csv")
glb_out_obs <- NULL # select from c(NULL : default to "new", "all", "new", "trn")
glb_out_vars_lst <- list()
# glb_id_var will be the first output column, by default
glb_out_vars_lst[["Probability1"]] <- 
    "%<d-% paste0(glb_rsp_var_out, glb_fin_mdl_id, '.prob')"
# glb_out_vars_lst[[glb_rsp_var_raw]] <- glb_rsp_var_raw
# glb_out_vars_lst[[paste0(head(unlist(strsplit(glb_rsp_var_out, "")), -1), collapse = "")]] <-

glbOutStackFnames <- NULL #: default
    # c("ebayipads_txt_assoc1_out_bid1_stack.csv") # manual stack
    # c("ebayipads_finmdl_bid1_out_nnet_1.csv") # universal stack
glb_out_pfx <- "NYTBlogs3_base_"
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

![](NYTBlogs3_base_files/figure-html/set_global_options-1.png) 

```r
glb_analytics_avl_objs <- NULL

glb_chunks_df <- myadd_chunk(NULL, "import.data")
```

```
##         label step_major step_minor label_minor   bgn end elapsed
## 1 import.data          1          0           0 9.037  NA      NA
```

## Step `1.0: import data`
#### chunk option: eval=<r condition>

```r
#glb_chunks_df <- myadd_chunk(NULL, "import.data")

glb_to_sav <- function() {
    sav_allobs_df <<- glbObsAll 
    sav_trnobs_df <<- glbObsTrn
    if (any(grepl("glbObsFit", ls(envir=globalenv()), fixed=TRUE)) &&
        !is.null(glbObsFit)) sav_fitobs_df <<- glbObsFit    
    if (any(grepl("glbObsOOB", ls(envir=globalenv()), fixed=TRUE)) &&
        !is.null(glbObsOOB)) sav_OOBobs_df <<- glbObsOOB    
    if (any(grepl("glbObsNew", ls(envir=globalenv()), fixed=TRUE)) &&
        !is.null(glbObsNew)) {
        #print("Attempting to save glbObsNew...")
        sav_newobs_df <<- glbObsNew    
    }

    if (any(grepl("glb_ctgry_df", ls(envir=globalenv()), fixed=TRUE)) &&
        !is.null(glb_ctgry_df)) sav_ctgry_df <<- glb_ctgry_df    

    if (!is.null(glb_models_lst )) sav_models_lst  <<- glb_models_lst
    if (!is.null(glb_models_df  )) sav_models_df   <<- glb_models_df

    if (any(grepl("glb_feats_df", ls(envir=globalenv()), fixed=TRUE)) &&
        !is.null(glb_feats_df)) sav_feats_df <<- glb_feats_df    
    if (any(grepl("glb_featsimp_df", ls(envir=globalenv()), fixed=TRUE)) &&
        !is.null(glb_featsimp_df)) sav_featsimp_df <<- glb_featsimp_df    
}

glbObsTrn <- myimport_data(url=glb_trnng_url, comment="glbObsTrn", 
                                force_header=TRUE)
```

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

```r
# glbObsTrn <- read.delim("data/hygiene.txt", header=TRUE, fill=TRUE, sep="\t",
#                             fileEncoding='iso-8859-1')
# glbObsTrn <- read.table("data/hygiene.dat.labels", col.names=c("dirty"),
#                             na.strings="[none]")
# glbObsTrn$review <- readLines("data/hygiene.dat", n =-1)
# comment(glbObsTrn) <- "glbObsTrn"                                

# glbObsTrn <- data.frame()
# for (symbol in c("Boeing", "CocaCola", "GE", "IBM", "ProcterGamble")) {
#     sym_trnobs_df <- 
#         myimport_data(url=gsub("IBM", symbol, glb_trnng_url), comment="glbObsTrn", 
#                                     force_header=TRUE)
#     sym_trnobs_df$Symbol <- symbol
#     glbObsTrn <- myrbind_df(glbObsTrn, sym_trnobs_df)
# }
                                
# glbObsTrn <- 
#     glbObsTrn %>% dplyr::filter(Year >= 1999)
                                
if (glb_is_separate_newobs_dataset) {
    glbObsNew <- myimport_data(url=glb_newdt_url, comment="glbObsNew", 
                                   force_header=TRUE)
    
    # To make plots / stats / checks easier in chunk:inspectORexplore.data
    glbObsAll <- myrbind_df(glbObsTrn, glbObsNew); 
    comment(glbObsAll) <- "glbObsAll"
} else {
    glbObsAll <- glbObsTrn; comment(glbObsAll) <- "glbObsAll"
    if (!glb_split_entity_newobs_datasets) {
        stop("Not implemented yet") 
        glbObsNew <- glbObsTrn[sample(1:nrow(glbObsTrn),
                                          max(2, nrow(glbObsTrn) / 1000)),]                    
    } else      if (glb_split_newdata_method == "condition") {
            glbObsNew <- do.call("subset", 
                list(glbObsTrn, parse(text=glb_split_newdata_condition)))
            glbObsTrn <- do.call("subset", 
                list(glbObsTrn, parse(text=paste0("!(", 
                                                      glb_split_newdata_condition,
                                                      ")"))))
        } else if (glb_split_newdata_method == "sample") {
                require(caTools)
                
                set.seed(glb_split_sample.seed)
                split <- sample.split(glbObsTrn[, glb_rsp_var_raw], 
                                      SplitRatio=(1-glb_split_newdata_size_ratio))
                glbObsNew <- glbObsTrn[!split, ] 
                glbObsTrn <- glbObsTrn[split ,]
        } else if (glb_split_newdata_method == "copy") {  
            glbObsTrn <- glbObsAll
            comment(glbObsTrn) <- "glbObsTrn"
            glbObsNew <- glbObsAll
            comment(glbObsNew) <- "glbObsNew"
        } else stop("glb_split_newdata_method should be %in% c('condition', 'sample', 'copy')")   

    comment(glbObsNew) <- "glbObsNew"
    myprint_df(glbObsNew)
    str(glbObsNew)

    if (glb_split_entity_newobs_datasets) {
        myprint_df(glbObsTrn)
        str(glbObsTrn)        
    }
}         
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

```r
if ((num_nas <- sum(is.na(glbObsTrn[, glb_rsp_var_raw]))) > 0)
    stop("glbObsTrn$", glb_rsp_var_raw, " contains NAs for ", num_nas, " obs")

if (nrow(glbObsTrn) == nrow(glbObsAll))
    warning("glbObsTrn same as glbObsAll")
if (nrow(glbObsNew) == nrow(glbObsAll))
    warning("glbObsNew same as glbObsAll")

if (length(glb_drop_vars) > 0) {
    warning("dropping vars: ", paste0(glb_drop_vars, collapse=", "))
    glbObsAll <- glbObsAll[, setdiff(names(glbObsAll), glb_drop_vars)]
    glbObsTrn <- glbObsTrn[, setdiff(names(glbObsTrn), glb_drop_vars)]    
    glbObsNew <- glbObsNew[, setdiff(names(glbObsNew), glb_drop_vars)]    
}

#stop(here"); sav_allobs_df <- glbObsAll # glbObsAll <- sav_allobs_df
# Combine trnent & newobs into glbObsAll for easier manipulation
glbObsTrn$.src <- "Train"; glbObsNew$.src <- "Test"; 
glbFeatsExclude <- union(glbFeatsExclude, ".src")
glbObsAll <- myrbind_df(glbObsTrn, glbObsNew)
comment(glbObsAll) <- "glbObsAll"

# Check for duplicates in glb_id_var
if (length(glb_id_var) == 0) {
    warning("using .rownames as identifiers for observations")
    glbObsAll$.rownames <- rownames(glbObsAll)
    glbObsTrn$.rownames <- rownames(subset(glbObsAll, .src == "Train"))
    glbObsNew$.rownames <- rownames(subset(glbObsAll, .src == "Test"))    
    glb_id_var <- ".rownames"
}
if (sum(duplicated(glbObsAll[, glb_id_var, FALSE])) > 0)
    stop(glb_id_var, " duplicated in glbObsAll")
glbFeatsExclude <- union(glbFeatsExclude, glb_id_var)

glbObsAll <- orderBy(reformulate(glb_id_var), glbObsAll)
glbObsTrn <- glbObsNew <- NULL

# For Tableau
if (!is.null(glbOutDataVizFname))
    write.csv(glbObsAll, glbOutDataVizFname, row.names=FALSE)

# - Merge glb_obs_stack_condition & glbObsDropCondition
# - Derive glb_obs_stack|drop_chk_vars from condition automatically
# - Implement glb_obs_stack_condition & glb_obs_stack_chk_vars options

dsp_partition_stats <- function(obs_df, vars=NULL) {
    
    lcl_vars <- NULL
    for (var in c(vars, glb_rsp_var_raw)) {
        if ((length(unique(obs_df[, var])) > 5) && is.numeric(obs_df[, var])) {
            cut_var <- paste0(var, ".cut.fctr")
            obs_df[, cut_var] <- cut(obs_df[, var], 3)
            lcl_vars <- union(lcl_vars, cut_var)
        } else lcl_vars <- union(lcl_vars, var)   
    }

    print("Partition stats:")
    print(mycreate_sqlxtab_df(obs_df, union(lcl_vars, ".src")))
    for (var in lcl_vars) {
        print(freq_df <- mycreate_sqlxtab_df(obs_df, union(var, ".src")))
        print(myplot_hbar(freq_df, ".src", ".n", colorcol_name=var))
    }
    print(mycreate_sqlxtab_df(obs_df, ".src"))
        
}

myget_symbols <- function(txt) {
    if (is.null(txt)) return(NULL)
    #print(getParseData(parse(text=txt, keep.source=TRUE)))
    return(unique(subset(getParseData(parse(text=txt, keep.source=TRUE)), 
                         token == "SYMBOL")$text))
}
# tokens <- unlist(strsplit(gsub("[[:punct:]|[:space:]]", " ", glbObsDropCondition), " "))
# tokens <- tokens[tokens != ""]
# glb_obs_drop_chk_vars <- c("biddable") # or NULL

dsp_partition_stats(obs_df=glbObsAll, vars=myget_symbols(glbObsDropCondition))
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

![](NYTBlogs3_base_files/figure-html/import.data-1.png) 

```
##    .src   .n
## 1 Train 6532
## 2  Test 1870
```

```r
if (!is.null(glbObsDropCondition)) {
    print(sprintf("Running glbObsDropCondition filter: %s", glbObsDropCondition))
    glbObsAll <- do.call("subset", 
                list(glbObsAll, parse(text=paste0("!(", glbObsDropCondition, ")"))))
    dsp_partition_stats(obs_df=glbObsAll, vars=myget_symbols(glbObsDropCondition))    
}

# Check for duplicates by all features
# Refactor to utilize glbSpecs
require(dplyr)
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
```

```r
require(lazyeval)
```

```
## Loading required package: lazyeval
```

```r
require(gdata)
```

```
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

```r
#print(names(glbObsAll))
dupObsAll <- glbObsAll[duplicated2(subset(glbObsAll, 
                                        select = -c(UniqueID, Popular, .src))), ]
dupObsAllIx <- glbObsAll %>% 
    dplyr::select_(.dots = setdiff(names(glbObsAll), 
                                   c(glb_id_var, glb_rsp_var_raw, ".src"))) %>%    
    #dplyr::select(-c(UniqueID, Popular, .src)) %>%
    duplicated2()
dupObsAll <- glbObsAll[dupObsAllIx, ] %>% 
    dplyr::arrange_(.dots = setdiff(names(glbObsAll), 
                                   c(glb_id_var, glb_rsp_var_raw, ".src")))
print(sprintf("Found %d duplicates by all features:", nrow(dupObsAll)))
```

```
## [1] "Found 0 duplicates by all features:"
```

```r
myprint_df(dupObsAll)
```

```
##  [1] NewsDesk       SectionName    SubsectionName Headline      
##  [5] Snippet        Abstract       WordCount      PubDate       
##  [9] Popular        UniqueID       .src          
## <0 rows> (or 0-length row.names)
```

```r
# print(dupObsAll[, c(glb_id_var, glb_rsp_var_raw, 
#                          "description", "startprice", "biddable")])
# write.csv(dupObsAll[, c("UniqueID"), FALSE], "ebayipads_dups.csv", row.names=FALSE)

if (nrow(dupObsAll) > 0) {
    dupobs_df <- tidyr::unite(dupObsAll, "allfeats", -c(sold, UniqueID, .src), sep="#")
    # dupobs_df <- dplyr::group_by(dupobs_df, allfeats)
    # dupobs_df <- dupobs_df[, "UniqueID", FALSE]
    # dupobs_df <- ungroup(dupobs_df)
    # 
    # dupobs_df$.rownames <- row.names(dupobs_df)
    grpobs_df <- data.frame(allfeats=unique(dupobs_df[, "allfeats"]))
    grpobs_df$.grpid <- row.names(grpobs_df)
    dupobs_df <- merge(dupobs_df, grpobs_df)
    
    # dupobs_tbl <- table(dupobs_df$.grpid)
    # print(max(dupobs_tbl))
    # print(dupobs_tbl[which.max(dupobs_tbl)])
    # print(dupobs_df[dupobs_df$.grpid == names(dupobs_tbl[which.max(dupobs_tbl)]), ])
    # print(dupobs_df[dupobs_df$.grpid == 106, ])
    # for (grpid in c(9, 17, 31, 36, 53))
    #     print(dupobs_df[dupobs_df$.grpid == grpid, ])
    dupgrps_df <- as.data.frame(table(dupobs_df$.grpid, dupobs_df$sold, useNA="ifany"))
    names(dupgrps_df)[c(1,2)] <- c(".grpid", "sold")
    dupgrps_df$.grpid <- as.numeric(as.character(dupgrps_df$.grpid))
    dupgrps_df <- tidyr::spread(dupgrps_df, sold, Freq)
    names(dupgrps_df)[-1] <- paste("sold", names(dupgrps_df)[-1], sep=".")
    dupgrps_df$.freq <- sapply(1:nrow(dupgrps_df), function(row) sum(dupgrps_df[row, -1]))
    myprint_df(orderBy(~-.freq, dupgrps_df))
    
    print("sold Conflicts:")
    print(subset(dupgrps_df, (sold.0 > 0) & (sold.1 > 0)))
    #dupobs_df[dupobs_df$.grpid == 4, ]
    glbObsAll <- merge(glbObsAll, dupobs_df[, c(glb_id_var, ".grpid")], 
                           by=glb_id_var, all.x=TRUE)
    if (nrow(subset(dupgrps_df, (sold.0 > 0) & (sold.1 > 0) & (sold.0 != sold.1))) > 0)
        stop("Duplicate conflicts are resolvable")
    #subset(glbObsAll, .grpid %in% c(25))
    #mydspObs(list(productline.contains="iPad 1", storage.contains="16", color.contains="Black", carrier.contains="None", cellular.contains="0", condition.contains="Used", startprice=80), cols=c("productline", "storage", "color", "carrier", "cellular", "condition", "startprice", "sold"))
    
    print("Test & Train Groups:")
    print(subset(dupgrps_df, (sold.NA > 0)))
    
    glbFeatsExclude <- c(".grpid", glbFeatsExclude)
}

if (!is.null(glbInpMerge)) {
    print("Running glbInpMerge specs...")
    obsMrg <- data.frame()
    for (fName in glbInpMerge$fnames) {
        print(sprintf("    Appending rows from %s...", fName))
        obsMrg <- rbind(obsMrg, read.csv(fName))
    }
    glbObsAll <- merge(glbObsAll, obsMrg, all.x = TRUE)
}

dsp_partition_stats(obs_df = glbObsAll,
                    vars = myget_symbols(glb_obs_repartition_train_condition))
```

```
## [1] "Partition stats:"
##   Popular  .src   .n
## 1       0 Train 5439
## 2      NA  Test 1870
## 3       1 Train 1093
##   Popular  .src   .n
## 1       0 Train 5439
## 2      NA  Test 1870
## 3       1 Train 1093
```

![](NYTBlogs3_base_files/figure-html/import.data-2.png) 

```
##    .src   .n
## 1 Train 6532
## 2  Test 1870
```

```r
if (!is.null(glb_obs_repartition_train_condition)) {
    print(sprintf("Running glb_obs_repartition_train_condition filter: %s",
                  glb_obs_repartition_train_condition))
#     glbObsAll <- mutate(glbObsAll, .src=ifelse(!is.na(sold) & (sold == 1),
#                             "Train", "Test"))
#     glbObsAll <- mutate_(glbObsAll, 
#                         .src=interp(ifelse(eval(parse(text="!is.na(sold) & (sold == 1)")),
#                                         "Train", "Test")))
#     glbObsAll <- within(glbObsAll, {
#         .src <- ifelse(eval(parse(text="!is.na(sold) & (sold == 1)")),
#                                         "Train", "Test")
#     })
#     glbObsAll <- within(glbObsAll, {
#         if(eval(parse(text="!is.na(sold) & (sold == 1)"))) .src <- "Train" else
#             .src <- "Test"
#     })
#     with(glbObsAll, {
#         src <- ifelse(eval(parse(text="!is.na(sold) & (sold == 1)")),
#                                         "Train", "Test")
#     })
#     glbObsAll$.src <- sapply(1:nrow(glbObsAll), function (row_ix) ifelse)
#     glbObsAll[parse(text=paste0("!(", glbObsDropCondition, ")")), ".src"] <- do.call("subset", 
#                 list(glbObsAll, ))
    
    glbObsTrn <- do.call("subset", list(glbObsAll, 
                        parse(text=paste0(" (", glb_obs_repartition_train_condition, ")"))))
    glbObsTrn$.src <- "Train"
    glbObsNew <- do.call("subset", list(glbObsAll, 
                        parse(text=paste0("!(", glb_obs_repartition_train_condition, ")"))))
    glbObsNew$.src <- "Test"
    glbObsAll <- rbind(glbObsTrn, glbObsNew)

    dsp_partition_stats(obs_df = glbObsAll,
                        vars = myget_symbols(glb_obs_repartition_train_condition))    
}

glb_chunks_df <- myadd_chunk(glb_chunks_df, "inspect.data", major.inc=TRUE)
```

```
##          label step_major step_minor label_minor    bgn    end elapsed
## 1  import.data          1          0           0  9.037 18.458   9.421
## 2 inspect.data          2          0           0 18.458     NA      NA
```

## Step `2.0: inspect data`

```r
#print(str(glbObsAll))
#View(glbObsAll)

dsp_class_dstrb <- function(var) {
    xtab_df <- mycreate_xtab_df(glbObsAll, c(".src", var))
    rownames(xtab_df) <- xtab_df$.src
    xtab_df <- subset(xtab_df, select=-.src)
    print(xtab_df)
    print(xtab_df / rowSums(xtab_df, na.rm=TRUE))    
}    

# Performed repeatedly in other chunks
glb_chk_data <- function(featsExclude = glbFeatsExclude, 
                         fctrMaxUniqVals = glbFctrMaxUniqVals) {
    # Histogram of predictor in glbObsTrn & glbObsNew
    print(myplot_histogram(glbObsAll, glb_rsp_var_raw) + facet_wrap(~ .src))
    
    if (glb_is_classification) 
        dsp_class_dstrb(var=ifelse(glb_rsp_var %in% names(glbObsAll), 
                                   glb_rsp_var, glb_rsp_var_raw))
    mycheck_problem_data(glbObsAll, featsExclude, fctrMaxUniqVals)
}
glb_chk_data()
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

```
## Warning: Removed 1870 rows containing non-finite values (stat_bin).
```

```
## Loading required package: reshape2
```

![](NYTBlogs3_base_files/figure-html/inspect.data-1.png) 

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

```r
# Create new features that help diagnostics
if (!is.null(glb_map_rsp_raw_to_var)) {
    glbObsAll[, glb_rsp_var] <- 
        glb_map_rsp_raw_to_var(glbObsAll[, glb_rsp_var_raw])
    mycheck_map_results(mapd_df=glbObsAll, 
                        from_col_name=glb_rsp_var_raw, to_col_name=glb_rsp_var)
        
    if (glb_is_classification) dsp_class_dstrb(glb_rsp_var)
}
```

```
##   Popular Popular.fctr   .n
## 1       0            N 5439
## 2      NA         <NA> 1870
## 3       1            Y 1093
```

```
## Warning: Removed 1 rows containing missing values (position_stack).
```

![](NYTBlogs3_base_files/figure-html/inspect.data-2.png) 

```
##       Popular.fctr.N Popular.fctr.Y Popular.fctr.NA
## Test              NA             NA            1870
## Train           5439           1093              NA
##       Popular.fctr.N Popular.fctr.Y Popular.fctr.NA
## Test              NA             NA               1
## Train      0.8326699      0.1673301              NA
```

```r
# check distribution of all numeric data
dsp_numeric_feats_dstrb <- function(feats_vctr) {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(ceiling(length(feats_vctr) / 2.0), 2)))
    pltIx <- 1
    for (feat in feats_vctr) {
        #print(sprintf("feat: %s", feat))
        if (glb_is_regression)
            gp <- myplot_scatter(df=glbObsAll, ycol_name=glb_rsp_var, xcol_name=feat,
                                 smooth=TRUE)
        if (glb_is_classification)
            #gp <- myplot_box(df=glbObsAll, ycol_names=feat, xcol_name=glb_rsp_var)
            gp <- myplot_violin(glbObsAll, ycol_names = feat, xcol_name = glb_rsp_var)
        if (inherits(glbObsAll[, feat], "factor"))
            gp <- gp + facet_wrap(reformulate(feat))
        print(gp + labs(title = feat), 
              vp = viewport(layout.pos.row = ceiling(pltIx / 2.0), 
                            layout.pos.col = ((pltIx - 1) %% 2) + 1))  
        
        pltIx <- pltIx + 1        
    }
}
# dsp_numeric_feats_dstrb(setdiff(names(glbObsAll), union(myfind_chr_cols_df(glbObsAll), c(glb_rsp_var_raw, glb_rsp_var)))[2:6])              # dsp_numeric_feats_dstrb(c("startprice", "sprice.root3", "sprice.predict.diff"))                                      
add_new_diag_feats <- function(obs_df, ref_df=glbObsAll) {
    require(plyr)
    
    set.seed(169)
    obs_df <- mutate(obs_df,
#         <col_name>.NA=is.na(<col_name>),

#         <col_name>.fctr=factor(<col_name>, 
#                     as.factor(union(obs_df$<col_name>, obs_twin_df$<col_name>))), 
#         <col_name>.fctr=relevel(factor(<col_name>, 
#                     as.factor(union(obs_df$<col_name>, obs_twin_df$<col_name>))),
#                                   "<ref_val>"), 
#         <col2_name>.fctr=relevel(factor(ifelse(<col1_name> == <val>, "<oth_val>", "<ref_val>")), 
#                               as.factor(c("R", "<ref_val>")),
#                               ref="<ref_val>"),

          # This doesn't work - use sapply instead
#         <col_name>.fctr_num=grep(<col_name>, levels(<col_name>.fctr)), 
#         
#         Date.my=as.Date(strptime(Date, "%m/%d/%y %H:%M")),
#         Year=year(Date.my),
#         Month=months(Date.my),
#         Weekday=weekdays(Date.my)

#         <col_name>=<table>[as.character(<col2_name>)],
#         <col_name>=as.numeric(<col2_name>),

#         <col_name> = trunc(<col2_name> / 100),

        .rnorm = rnorm(n=nrow(obs_df))
                        )

    # If levels of a factor are different across obs_df & glbObsNew; predict.glm fails  
    # Transformations not handled by mutate
#     obs_df$<col_name>.fctr.num <- sapply(1:nrow(obs_df), 
#         function(row_ix) grep(obs_df[row_ix, "<col_name>"],
#                               levels(obs_df[row_ix, "<col_name>.fctr"])))
    
    #print(summary(obs_df))
    #print(sapply(names(obs_df), function(col) sum(is.na(obs_df[, col]))))
    return(obs_df)
}
glbObsAll <- add_new_diag_feats(glbObsAll)
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
```

```r
require(dplyr)

#stop(here"); sav_allobs_df <- glbObsAll # glbObsAll <- sav_allobs_df
# Merge some <descriptor>
# glbObsAll$<descriptor>.my <- glbObsAll$<descriptor>
# glbObsAll[grepl("\\bAIRPORT\\b", glbObsAll$<descriptor>.my),
#               "<descriptor>.my"] <- "AIRPORT"

# Check distributions of newly transformed / extracted vars
#   Enhancement: remove vars that were displayed ealier
dsp_numeric_feats_dstrb(feats_vctr=setdiff(names(glbObsAll), 
        c(myfind_chr_cols_df(glbObsAll), glb_rsp_var_raw, glb_rsp_var, 
          glbFeatsExclude)))
```

![](NYTBlogs3_base_files/figure-html/inspect.data-3.png) 

```r
#   Convert factors to dummy variables
#   Build splines   require(splines); bsBasis <- bs(training$age, df=3)

#pairs(subset(glbObsTrn, select=-c(col_symbol)))
# Check for glbObsNew & glbObsTrn features range mismatches

# Other diagnostics:
# print(subset(glbObsTrn, <col1_name> == max(glbObsTrn$<col1_name>, na.rm=TRUE) & 
#                         <col2_name> <= mean(glbObsTrn$<col1_name>, na.rm=TRUE)))

# print(glbObsTrn[which.max(glbObsTrn$<col_name>),])

# print(<col_name>_freq_glbObsTrn <- mycreate_tbl_df(glbObsTrn, "<col_name>"))
# print(which.min(table(glbObsTrn$<col_name>)))
# print(which.max(table(glbObsTrn$<col_name>)))
# print(which.max(table(glbObsTrn$<col1_name>, glbObsTrn$<col2_name>)[, 2]))
# print(table(glbObsTrn$<col1_name>, glbObsTrn$<col2_name>))
# print(table(is.na(glbObsTrn$<col1_name>), glbObsTrn$<col2_name>))
# print(table(sign(glbObsTrn$<col1_name>), glbObsTrn$<col2_name>))
# print(mycreate_xtab_df(glbObsTrn, <col1_name>))
# print(mycreate_xtab_df(glbObsTrn, c(<col1_name>, <col2_name>)))
# print(<col1_name>_<col2_name>_xtab_glbObsTrn <- 
#   mycreate_xtab_df(glbObsTrn, c("<col1_name>", "<col2_name>")))
# <col1_name>_<col2_name>_xtab_glbObsTrn[is.na(<col1_name>_<col2_name>_xtab_glbObsTrn)] <- 0
# print(<col1_name>_<col2_name>_xtab_glbObsTrn <- 
#   mutate(<col1_name>_<col2_name>_xtab_glbObsTrn, 
#             <col3_name>=(<col1_name> * 1.0) / (<col1_name> + <col2_name>))) 
# print(mycreate_sqlxtab_df(glbObsAll, c("<col1_name>", "<col2_name>")))

# print(<col2_name>_min_entity_arr <- 
#    sort(tapply(glbObsTrn$<col1_name>, glbObsTrn$<col2_name>, min, na.rm=TRUE)))
# print(<col1_name>_na_by_<col2_name>_arr <- 
#    sort(tapply(glbObsTrn$<col1_name>.NA, glbObsTrn$<col2_name>, mean, na.rm=TRUE)))

# Other plots:
# print(myplot_box(df=glbObsTrn, ycol_names="<col1_name>"))
# print(myplot_box(df=glbObsTrn, ycol_names="<col1_name>", xcol_name="<col2_name>"))
# print(myplot_line(subset(glbObsTrn, Symbol %in% c("CocaCola", "ProcterGamble")), 
#                   "Date.POSIX", "StockPrice", facet_row_colnames="Symbol") + 
#     geom_vline(xintercept=as.numeric(as.POSIXlt("2003-03-01"))) +
#     geom_vline(xintercept=as.numeric(as.POSIXlt("1983-01-01")))        
#         )
# print(myplot_line(subset(glbObsTrn, Date.POSIX > as.POSIXct("2004-01-01")), 
#                   "Date.POSIX", "StockPrice") +
#     geom_line(aes(color=Symbol)) + 
#     coord_cartesian(xlim=c(as.POSIXct("1990-01-01"),
#                            as.POSIXct("2000-01-01"))) +     
#     coord_cartesian(ylim=c(0, 250)) +     
#     geom_vline(xintercept=as.numeric(as.POSIXlt("1997-09-01"))) +
#     geom_vline(xintercept=as.numeric(as.POSIXlt("1997-11-01")))        
#         )
# print(myplot_scatter(glbObsAll, "<col1_name>", "<col2_name>", smooth=TRUE))
# print(myplot_scatter(glbObsAll, "<col1_name>", "<col2_name>", colorcol_name="<Pred.fctr>") + 
#         geom_point(data=subset(glbObsAll, <condition>), 
#                     mapping=aes(x=<x_var>, y=<y_var>), color="red", shape=4, size=5) +
#         geom_vline(xintercept=84))

glb_chunks_df <- myadd_chunk(glb_chunks_df, "scrub.data", major.inc=FALSE)
```

```
##          label step_major step_minor label_minor    bgn   end elapsed
## 2 inspect.data          2          0           0 18.458 21.34   2.882
## 3   scrub.data          2          1           1 21.340    NA      NA
```

### Step `2.1: scrub data`

```r
mycheck_problem_data(glbObsAll, featsExclude = glbFeatsExclude, 
                     fctrMaxUniqVals = glbFctrMaxUniqVals)
```

```
## [1] "numeric data missing in glbObsAll: "
##      Popular Popular.fctr 
##         1870         1870 
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

```r
findOffendingCharacter <- function(x, maxStringLength=256){  
  print(x)
  for (c in 1:maxStringLength){
    offendingChar <- substr(x,c,c)
    #print(offendingChar) #uncomment if you want the indiv characters printed
    #the next character is the offending multibyte Character
  }    
}
# string_vector <- c("test", "Se\x96ora", "works fine")
# lapply(string_vector, findOffendingCharacter)
# lapply(glbObsAll$description[29], findOffendingCharacter)

dsp_hdlxtab <- function(str) 
    print(mycreate_sqlxtab_df(glbObsAll[sel_obs(Headline.contains=str), ],
                           c("Headline.pfx", "Headline", glb_rsp_var)))
#dsp_hdlxtab("(1914)|(1939)")

dsp_catxtab <- function(str) 
    print(mycreate_sqlxtab_df(glbObsAll[sel_obs(Headline.contains=str), ],
        c("Headline.pfx", "NewsDesk", "SectionName", "SubsectionName", glb_rsp_var)))
# dsp_catxtab("1914)|(1939)")
# dsp_catxtab("19(14|39|64):")
# dsp_catxtab("19..:")

# Merge some categories
# glbObsAll$myCategory <-
#     plyr::revalue(glbObsAll$myCategory, c(      
#         "#Business Day#Dealbook"            = "Business#Business Day#Dealbook",
#         "#Business Day#Small Business"      = "Business#Business Day#Small Business",
#         "dummy" = "dummy"
#     ))

# ctgry_xtab_df <- orderBy(reformulate(c("-", ".n")),
#                           mycreate_sqlxtab_df(glbObsAll,
#     c("myCategory", "NewsDesk", "SectionName", "SubsectionName", glb_rsp_var)))
# myprint_df(ctgry_xtab_df)
# write.table(ctgry_xtab_df, paste0(glb_out_pfx, "ctgry_xtab.csv"), 
#             row.names=FALSE)

# ctgry_cast_df <- orderBy(~ -Y -NA, dcast(ctgry_xtab_df, 
#                        myCategory + NewsDesk + SectionName + SubsectionName ~ 
#                            Popular.fctr, sum, value.var=".n"))
# myprint_df(ctgry_cast_df)
# write.table(ctgry_cast_df, paste0(glb_out_pfx, "ctgry_cast.csv"), 
#             row.names=FALSE)

# print(ctgry_sum_tbl <- table(glbObsAll$myCategory, glbObsAll[, glb_rsp_var], 
#                              useNA="ifany"))

dsp_chisq.test <- function(...) {
    sel_df <- glbObsAll[sel_obs(...) & 
                            !is.na(glbObsAll$Popular), ]
    sel_df$.marker <- 1
    ref_df <- glbObsAll[!is.na(glbObsAll$Popular), ]
    mrg_df <- merge(ref_df[, c(glb_id_var, "Popular")],
                    sel_df[, c(glb_id_var, ".marker")], all.x=TRUE)
    mrg_df[is.na(mrg_df)] <- 0
    print(mrg_tbl <- table(mrg_df$.marker, mrg_df$Popular))
    print("Rows:Selected; Cols:Popular")
    #print(mrg_tbl)
    print(chisq.test(mrg_tbl))
}
# dsp_chisq.test(Headline.contains="[Ee]bola")
# dsp_chisq.test(Snippet.contains="[Ee]bola")
# dsp_chisq.test(Abstract.contains="[Ee]bola")

# print(mycreate_sqlxtab_df(glbObsAll[sel_obs(Headline.contains="[Ee]bola"), ], 
#                           c(glb_rsp_var, "NewsDesk", "SectionName", "SubsectionName")))

# print(table(glbObsAll$NewsDesk, glbObsAll$SectionName))
# print(table(glbObsAll$SectionName, glbObsAll$SubsectionName))
# print(table(glbObsAll$NewsDesk, glbObsAll$SectionName, glbObsAll$SubsectionName))

# glbObsAll$myCategory.fctr <- as.factor(glbObsAll$myCategory)
```

### Step `2.1: scrub data`

```r
glb_chunks_df <- myadd_chunk(glb_chunks_df, "transform.data", major.inc=FALSE)
```

```
##            label step_major step_minor label_minor    bgn   end elapsed
## 3     scrub.data          2          1           1 21.340 22.44     1.1
## 4 transform.data          2          2           2 22.441    NA      NA
```

```r
### Mapping dictionary
#sav_allobs_df <- glbObsAll; glbObsAll <- sav_allobs_df
if (!is.null(glb_map_vars)) {
    for (feat in glb_map_vars) {
        map_df <- myimport_data(url=glb_map_urls[[feat]], 
                                            comment="map_df", 
                                           print_diagn=TRUE)
        glbObsAll <- mymap_codes(glbObsAll, feat, names(map_df)[2], 
                                     map_df, map_join_col_name=names(map_df)[1], 
                                     map_tgt_col_name=names(map_df)[2])
    }
    glbFeatsExclude <- union(glbFeatsExclude, glb_map_vars)
}

### Forced Assignments
#stop(here"); sav_allobs_df <- glbObsAll; glbObsAll <- sav_allobs_df
for (feat in glb_assign_vars) {
    new_feat <- paste0(feat, ".my")
    print(sprintf("Forced Assignments for: %s -> %s...", feat, new_feat))
    glbObsAll[, new_feat] <- glbObsAll[, feat]
    
    pairs <- glb_assign_pairs_lst[[feat]]
    for (pair_ix in 1:length(pairs$from)) {
        if (is.na(pairs$from[pair_ix]))
            nobs <- nrow(filter(glbObsAll, 
                                is.na(eval(parse(text=feat),
                                            envir=glbObsAll)))) else
            nobs <- sum(glbObsAll[, feat] == pairs$from[pair_ix])
        #nobs <- nrow(filter(glbObsAll, is.na(Married.fctr)))    ; print(nobs)
        
        if ((is.na(pairs$from[pair_ix])) && (is.na(pairs$to[pair_ix])))
            stop("what are you trying to do ???")
        if (is.na(pairs$from[pair_ix]))
            glbObsAll[is.na(glbObsAll[, feat]), new_feat] <- 
                pairs$to[pair_ix] else
            glbObsAll[glbObsAll[, feat] == pairs$from[pair_ix], new_feat] <- 
                pairs$to[pair_ix]
                    
        print(sprintf("    %s -> %s for %s obs", 
                      pairs$from[pair_ix], pairs$to[pair_ix], format(nobs, big.mark=",")))
    }

    glbFeatsExclude <- union(glbFeatsExclude, glb_assign_vars)
}

### Derivations using mapping functions
#stop(here"); sav_allobs_df <- glbObsAll; glbObsAll <- sav_allobs_df
for (new_feat in glb_derive_vars) {
    print(sprintf("Creating new feature: %s...", new_feat))
    args_lst <- NULL 
    for (arg in glbFeatsDerive[[new_feat]]$args) 
        args_lst[[arg]] <- glbObsAll[, arg]
    glbObsAll[, new_feat] <- do.call(glbFeatsDerive[[new_feat]]$mapfn, args_lst)
}
```

```
## [1] "Creating new feature: NDSSName.my..."
## [1] "Creating new feature: WordCount.log1p..."
## [1] "Creating new feature: WordCount.root2..."
## [1] "Creating new feature: WordCount.nexp..."
```

```r
#stop(here")
#hex_vctr <- c("\n", "\211", "\235", "\317", "\333")
hex_regex <- paste0(c("\n", "\211", "\235", "\317", "\333"), collapse="|")
for (obs_id in c(10029, 10948, 10136, 10178, 11514, 11904, 12157, 12210, 12659)) {
#     tmp_str <- unlist(strsplit(glbObsAll[row_pos, "descr.my"], ""))
#     glbObsAll[row_pos, "descr.my"] <- paste0(tmp_str[!tmp_str %in% hex_vctr],
#                                                          collapse="")
    row_pos <- which(glbObsAll$UniqueID == obs_id)
#     glbObsAll[row_pos, "descr.my"] <- 
#         gsub(hex_regex, " ", glbObsAll[row_pos, "descr.my"])
}
```

## Step `2.2: transform data`

```r
#```{r extract_features, cache=FALSE, eval=!is.null(glbFeatsText)}
glb_chunks_df <- myadd_chunk(glb_chunks_df, "extract.features", major.inc=TRUE)
```

```
##              label step_major step_minor label_minor    bgn    end elapsed
## 4   transform.data          2          2           2 22.441 22.531   0.091
## 5 extract.features          3          0           0 22.532     NA      NA
```

```r
extract.features_chunk_df <- myadd_chunk(NULL, "extract.features_bgn")
```

```
##                  label step_major step_minor label_minor    bgn end
## 1 extract.features_bgn          1          0           0 22.539  NA
##   elapsed
## 1      NA
```

```r
# Create new features that help prediction
# <col_name>.lag.2 <- lag(zoo(glbObsTrn$<col_name>), -2, na.pad=TRUE)
# glbObsTrn[, "<col_name>.lag.2"] <- coredata(<col_name>.lag.2)
# <col_name>.lag.2 <- lag(zoo(glbObsNew$<col_name>), -2, na.pad=TRUE)
# glbObsNew[, "<col_name>.lag.2"] <- coredata(<col_name>.lag.2)
# 
# glbObsNew[1, "<col_name>.lag.2"] <- glbObsTrn[nrow(glbObsTrn) - 1, 
#                                                    "<col_name>"]
# glbObsNew[2, "<col_name>.lag.2"] <- glbObsTrn[nrow(glbObsTrn), 
#                                                    "<col_name>"]
                                                   
# glbObsAll <- mutate(glbObsAll,
#     A.P.http=ifelse(grepl("http",Added,fixed=TRUE), 1, 0)
#                     )
# 
# glbObsTrn <- mutate(glbObsTrn,
#                     )
# 
# glbObsNew <- mutate(glbObsNew,
#                     )

#   Convert dates to numbers 
#       typically, dates come in as chars; 
#           so this must be done before converting chars to factors

#stop(here"); sav_allobs_df <- glbObsAll #; glbObsAll <- sav_allobs_df
if (!is.null(glb_date_vars)) {
    glbObsAll <- cbind(glbObsAll, 
        myextract_dates_df(df=glbObsAll, vars=glb_date_vars, 
                           id_vars=glb_id_var, rsp_var=glb_rsp_var))
    for (sfx in c("", ".POSIX"))
        glbFeatsExclude <- 
            union(glbFeatsExclude, 
                    paste(glb_date_vars, sfx, sep=""))

    for (feat in glb_date_vars) {
        glbObsAll <- orderBy(reformulate(paste0(feat, ".POSIX")), glbObsAll)
#         print(myplot_scatter(glbObsAll, xcol_name=paste0(feat, ".POSIX"),
#                              ycol_name=glb_rsp_var, colorcol_name=glb_rsp_var))
        print(myplot_scatter(glbObsAll[glbObsAll[, paste0(feat, ".POSIX")] >=
                                               strptime("2012-12-01", "%Y-%m-%d"), ], 
                             xcol_name=paste0(feat, ".POSIX"),
                             ycol_name=glb_rsp_var, colorcol_name=paste0(feat, ".wkend")))

        # Create features that measure the gap between previous timestamp in the data
        require(zoo)
        z <- zoo(as.numeric(as.POSIXlt(glbObsAll[, paste0(feat, ".POSIX")])))
        glbObsAll[, paste0(feat, ".zoo")] <- z
        print(head(glbObsAll[, c(glb_id_var, feat, paste0(feat, ".zoo"))]))
        print(myplot_scatter(glbObsAll[glbObsAll[,  paste0(feat, ".POSIX")] >
                                            strptime("2012-10-01", "%Y-%m-%d"), ], 
                            xcol_name=paste0(feat, ".zoo"), ycol_name=glb_rsp_var,
                            colorcol_name=glb_rsp_var))
        b <- zoo(, seq(nrow(glbObsAll)))
        
        last1 <- as.numeric(merge(z-lag(z, -1), b, all=TRUE)); last1[is.na(last1)] <- 0
        glbObsAll[, paste0(feat, ".last1.log")] <- log(1 + last1)
        print(gp <- myplot_box(df=glbObsAll[glbObsAll[, 
                                                    paste0(feat, ".last1.log")] > 0, ], 
                               ycol_names=paste0(feat, ".last1.log"), 
                               xcol_name=glb_rsp_var))
        
        last2 <- as.numeric(merge(z-lag(z, -2), b, all=TRUE)); last2[is.na(last2)] <- 0
        glbObsAll[, paste0(feat, ".last2.log")] <- log(1 + last2)
        print(gp <- myplot_box(df=glbObsAll[glbObsAll[, 
                                                    paste0(feat, ".last2.log")] > 0, ], 
                               ycol_names=paste0(feat, ".last2.log"), 
                               xcol_name=glb_rsp_var))
        
        last10 <- as.numeric(merge(z-lag(z, -10), b, all=TRUE)); last10[is.na(last10)] <- 0
        glbObsAll[, paste0(feat, ".last10.log")] <- log(1 + last10)
        print(gp <- myplot_box(df=glbObsAll[glbObsAll[, 
                                                    paste0(feat, ".last10.log")] > 0, ], 
                               ycol_names=paste0(feat, ".last10.log"), 
                               xcol_name=glb_rsp_var))
        
        last100 <- as.numeric(merge(z-lag(z, -100), b, all=TRUE)); last100[is.na(last100)] <- 0
        glbObsAll[, paste0(feat, ".last100.log")] <- log(1 + last100)
        print(gp <- myplot_box(df=glbObsAll[glbObsAll[, 
                                                    paste0(feat, ".last100.log")] > 0, ], 
                               ycol_names=paste0(feat, ".last100.log"), 
                               xcol_name=glb_rsp_var))
        
        glbObsAll <- orderBy(reformulate(glb_id_var), glbObsAll)
        glbFeatsExclude <- union(glbFeatsExclude, 
                                                c(paste0(feat, ".zoo")))
        # all2$last3 = as.numeric(merge(z-lag(z, -3), b, all = TRUE))
        # all2$last5 = as.numeric(merge(z-lag(z, -5), b, all = TRUE))
        # all2$last10 = as.numeric(merge(z-lag(z, -10), b, all = TRUE))
        # all2$last20 = as.numeric(merge(z-lag(z, -20), b, all = TRUE))
        # all2$last50 = as.numeric(merge(z-lag(z, -50), b, all = TRUE))
        # 
        # 
        # # order table
        # all2 = all2[order(all2$id),]
        # 
        # ## fill in NAs
        # # count averages
        # na.avg = all2 %>% group_by(weekend, hour) %>% dplyr::summarise(
        #     last1=mean(last1, na.rm=TRUE),
        #     last3=mean(last3, na.rm=TRUE),
        #     last5=mean(last5, na.rm=TRUE),
        #     last10=mean(last10, na.rm=TRUE),
        #     last20=mean(last20, na.rm=TRUE),
        #     last50=mean(last50, na.rm=TRUE)
        # )
        # 
        # # fill in averages
        # na.merge = merge(all2, na.avg, by=c("weekend","hour"))
        # na.merge = na.merge[order(na.merge$id),]
        # for(i in c("last1", "last3", "last5", "last10", "last20", "last50")) {
        #     y = paste0(i, ".y")
        #     idx = is.na(all2[[i]])
        #     all2[idx,][[i]] <- na.merge[idx,][[y]]
        # }
        # rm(na.avg, na.merge, b, i, idx, n, pd, sec, sh, y, z)
    }
}
#rm(last1, last10, last100)

#   Create factors of string variables
extract.features_chunk_df <- myadd_chunk(extract.features_chunk_df, 
            paste0("extract.features_", "factorize.str.vars"), major.inc=TRUE)
```

```
##                                 label step_major step_minor label_minor
## 1                extract.features_bgn          1          0           0
## 2 extract.features_factorize.str.vars          2          0           0
##      bgn    end elapsed
## 1 22.539 22.552   0.013
## 2 22.552     NA      NA
```

```r
#stop(here"); sav_allobs_df <- glbObsAll; #glbObsAll <- sav_allobs_df
print(str_vars <- myfind_chr_cols_df(glbObsAll))
```

```
##         NewsDesk      SectionName   SubsectionName         Headline 
##       "NewsDesk"    "SectionName" "SubsectionName"       "Headline" 
##          Snippet         Abstract          PubDate             .src 
##        "Snippet"       "Abstract"        "PubDate"           ".src" 
##      NDSSName.my 
##    "NDSSName.my"
```

```r
if (length(str_vars <- setdiff(str_vars, 
                               c(glbFeatsExclude, glbFeatsText))) > 0) {
    for (var in str_vars) {
        warning("Creating factors of string variable: ", var, 
                ": # of unique values: ", length(unique(glbObsAll[, var])))
        glbObsAll[, paste0(var, ".fctr")] <- 
            relevel(factor(glbObsAll[, var]),
                    names(which.max(table(glbObsAll[, var], useNA = "ifany"))))
    }
    glbFeatsExclude <- union(glbFeatsExclude, str_vars)
}
```

```
## Warning: Creating factors of string variable: NDSSName.my: # of unique
## values: 21
```

```r
if (!is.null(glbFeatsText)) {
    require(foreach)
    require(gsubfn)
    require(stringr)
    require(tm)
    
    extract.features_chunk_df <- myadd_chunk(extract.features_chunk_df, 
            paste0("extract.features_", "process.text"), major.inc=TRUE)
    
    chk_pattern_freq <- function(rex_str, ignore.case=TRUE) {
        match_mtrx <- str_extract_all(txt_vctr, regex(rex_str, ignore_case=ignore.case), 
                                      simplify=TRUE)
        match_df <- as.data.frame(match_mtrx[match_mtrx != ""])
        names(match_df) <- "pattern"
        return(mycreate_sqlxtab_df(match_df, "pattern"))        
    }

#     match_lst <- gregexpr("\\bok(?!ay)", txt_vctr[746], ignore.case = FALSE, perl=TRUE); print(match_lst)
    dsp_pattern <- function(rex_str, ignore.case=TRUE, print.all=TRUE) {
        match_lst <- gregexpr(rex_str, txt_vctr, ignore.case = ignore.case, perl=TRUE)
        match_lst <- regmatches(txt_vctr, match_lst)
        match_df <- data.frame(matches=sapply(match_lst, 
                                              function (elems) paste(elems, collapse="#")))
        match_df <- subset(match_df, matches != "")
        if (print.all)
            print(match_df)
        return(match_df)
    }
    
    dsp_matches <- function(rex_str, ix) {
        print(match_pos <- gregexpr(rex_str, txt_vctr[ix], perl=TRUE))
        print(str_sub(txt_vctr[ix], (match_pos[[1]] / 100) *  99 +   0, 
                                    (match_pos[[1]] / 100) * 100 + 100))        
    }

    myapply_gsub <- function(...) {
        if ((length_lst <- length(names(gsub_map_lst))) == 0)
            return(txt_vctr)
        for (ptn_ix in 1:length_lst) {
            if ((ptn_ix %% 10) == 0)
                print(sprintf("running gsub for %02d (of %02d): #%s#...", ptn_ix, 
                                length(names(gsub_map_lst)), names(gsub_map_lst)[ptn_ix]))
            txt_vctr <- gsub(names(gsub_map_lst)[ptn_ix], gsub_map_lst[[ptn_ix]], 
                               txt_vctr, ...)
        }
        return(txt_vctr)
    }    

    myapply_txtmap <- function(txt_vctr, ...) {
        nrows <- nrow(glb_txt_map_df)
        for (ptn_ix in 1:nrows) {
            if ((ptn_ix %% 10) == 0)
                print(sprintf("running gsub for %02d (of %02d): #%s#...", ptn_ix, 
                                nrows, glb_txt_map_df[ptn_ix, "rex_str"]))
            txt_vctr <- gsub(glb_txt_map_df[ptn_ix, "rex_str"], 
                             glb_txt_map_df[ptn_ix, "rpl_str"], 
                               txt_vctr, ...)
        }
        return(txt_vctr)
        #print(txt_vctr <- glbObsAll[glbObsAll$UniqueID == 11329, "descr.my"])
        #strsplit(txt_vctr, "")[[1]][1]
        #ptn_ix <- 2; glb_txt_map_df[ptn_ix, ]
        #gsub(glb_txt_map_df[ptn_ix, "rex_str"], glb_txt_map_df[ptn_ix, "rpl_str"], txt_vctr)
        #print(match_lst <- gregexpr(glb_txt_map_df[ptn_ix, "rex_str"], txt_vctr))
        #strsplit(glb_txt_map_df[ptn_ix, "rex_str"], "")[[1]]
    }    

    chk.equal <- function(bgn, end) {
        print(all.equal(sav_txt_lst[["Headline"]][bgn:end], 
                        glb_txt_chr_lst[["Headline"]][bgn:end]))
    }    
    dsp.equal <- function(bgn, end) {
        print(sav_txt_lst[["Headline"]][bgn:end])
        print(glb_txt_chr_lst[["Headline"]][bgn:end])
    }    
#sav_txt_lst <- glb_txt_chr_lst; all.equal(sav_txt_lst, glb_txt_chr_lst)
#all.equal(sav_txt_lst[["Headline"]][1:4200], glb_txt_chr_lst[["Headline"]][1:4200])
#chk.equal( 1, 100)
#dsp.equal(86, 90)
    
    txt_map_filename <- paste0(glb_txt_munge_filenames_pfx, "map.csv")
    if (!file.exists(txt_map_filename))
        stop(txt_map_filename, " not found!")
    glb_txt_map_df <- read.csv(txt_map_filename, comment.char="#", strip.white=TRUE)
    glb_txt_chr_lst <- list(); 
    print(sprintf("Building glb_txt_chr_lst..."))
    glb_txt_chr_lst <- foreach(txt_var = glbFeatsText) %dopar% {   
#     for (txt_var in glbFeatsText) {
        txt_vctr <- glbObsAll[, txt_var]
        names(txt_vctr) <- glbObsAll[, glb_id_var]
        
        # myapply_txtmap shd be created as a tm_map::content_transformer ?
        #print(glb_txt_map_df)
        #txt_var=glbFeatsText[3]; txt_vctr <- glb_txt_chr_lst[[txt_var]]
        #print(rex_str <- glb_txt_map_df[3, "rex_str"])
        #print(rex_str <- glb_txt_map_df[glb_txt_map_df$rex_str == "\\bWall St\\.", "rex_str"])
        #print(rex_str <- glb_txt_map_df[grepl("du Pont", glb_txt_map_df$rex_str), "rex_str"])        
        #print(rex_str <- glb_txt_map_df[glb_txt_map_df$rpl_str == "versus", "rex_str"])             
        #print(tmp_vctr <- grep(rex_str, txt_vctr, value=TRUE, ignore.case=FALSE))
        #ret_lst <- regexec(rex_str, txt_vctr, ignore.case=FALSE); ret_lst <- regmatches(txt_vctr, ret_lst); ret_vctr <- sapply(1:length(ret_lst), function(pos_ix) ifelse(length(ret_lst[[pos_ix]]) > 0, ret_lst[[pos_ix]], "")); print(ret_vctr <- ret_vctr[ret_vctr != ""])
        #gsub(rex_str, glb_txt_map_df[glb_txt_map_df$rex_str == rex_str, "rpl_str"], tmp_vctr, ignore.case=FALSE)
        #grep("Hong Hong", txt_vctr, value=TRUE)
    
        txt_vctr <- myapply_txtmap(txt_vctr, ignore.case=FALSE)    
    }
    names(glb_txt_chr_lst) <- glbFeatsText

    for (txt_var in glbFeatsText) {
        print(sprintf("Remaining OK in %s:", txt_var))
        txt_vctr <- glb_txt_chr_lst[[txt_var]]
        
        print(chk_pattern_freq(rex_str <- "(?<!(BO|HO|LO))OK(?!(E\\!|ED|IE|IN|S ))",
                               ignore.case=FALSE))
        match_df <- dsp_pattern(rex_str, ignore.case=FALSE, print.all=FALSE)
        for (row in row.names(match_df))
            dsp_matches(rex_str, ix=as.numeric(row))

        print(chk_pattern_freq(rex_str <- "Ok(?!(a\\.|ay|in|ra|um))", ignore.case=FALSE))
        match_df <- dsp_pattern(rex_str, ignore.case=FALSE, print.all=FALSE)
        for (row in row.names(match_df))
            dsp_matches(rex_str, ix=as.numeric(row))

        print(chk_pattern_freq(rex_str <- "(?<!( b| B| c| C| g| G| j| M| p| P| w| W| r| Z|\\(b|ar|bo|Bo|co|Co|Ew|gk|go|ho|ig|jo|kb|ke|Ke|ki|lo|Lo|mo|mt|no|No|po|ra|ro|sm|Sm|Sp|to|To))ok(?!(ay|bo|e |e\\)|e,|e\\.|eb|ed|el|en|er|es|ey|i |ie|in|it|ka|ke|ki|ly|on|oy|ra|st|u |uc|uy|yl|yo))",
                               ignore.case=FALSE))
        match_df <- dsp_pattern(rex_str, ignore.case=FALSE, print.all=FALSE)
        for (row in row.names(match_df))
            dsp_matches(rex_str, ix=as.numeric(row))
    }    
    # txt_vctr <- glb_txt_chr_lst[[glbFeatsText[1]]]
    # print(chk_pattern_freq(rex_str <- "(?<!( b| c| C| p|\\(b|bo|co|lo|Lo|Sp|to|To))ok(?!(ay|e |e\\)|e,|e\\.|ed|el|en|es|ey|ie|in|on|ra))", ignore.case=FALSE))
    # print(chk_pattern_freq(rex_str <- "ok(?!(ay|el|on|ra))", ignore.case=FALSE))
    # dsp_pattern(rex_str, ignore.case=FALSE, print.all=FALSE)
    # dsp_matches(rex_str, ix=8)
    # substr(txt_vctr[86], 5613, 5620)
    # substr(glbObsAll[301, "review"], 550, 650)

#stop(here"); sav_txt_lst <- glb_txt_chr_lst    
    for (txt_var in glbFeatsText) {
        print(sprintf("Remaining Acronyms in %s:", txt_var))
        txt_vctr <- glb_txt_chr_lst[[txt_var]]
        
        print(chk_pattern_freq(rex_str <- "([[:upper:]]\\.( *)){2,}", ignore.case=FALSE))
        
        # Check for names
        print(subset(chk_pattern_freq(rex_str <- "(([[:upper:]]+)\\.( *)){1}",
                                      ignore.case=FALSE),
                     .n > 1))
        # dsp_pattern(rex_str="(OK\\.( *)){1}", ignore.case=FALSE)
        # dsp_matches(rex_str="(OK\\.( *)){1}", ix=557)
        #dsp_matches(rex_str="\\bR\\.I\\.P(\\.*)(\\B)", ix=461)
        #dsp_matches(rex_str="\\bR\\.I\\.P(\\.*)", ix=461)        
        #print(str_sub(txt_vctr[676], 10100, 10200))
        #print(str_sub(txt_vctr[74], 1, -1))        
    }

    for (txt_var in glbFeatsText) {
        re_str <- "\\b(Fort|Ft\\.|Hong|Las|Los|New|Puerto|Saint|San|St\\.)( |-)(\\w)+"
        print(sprintf("Remaining #%s# terms in %s: ", re_str, txt_var))
        txt_vctr <- glb_txt_chr_lst[[txt_var]]        
        print(orderBy(~ -.n +pattern, subset(chk_pattern_freq(re_str, ignore.case=FALSE), 
                                             grepl("( |-)[[:upper:]]", pattern))))
        print("    consider cleaning if relevant to problem domain; geography name; .n > 1")
        #grep("New G", txt_vctr, value=TRUE, ignore.case=FALSE)
        #grep("St\\. Wins", txt_vctr, value=TRUE, ignore.case=FALSE)
    }        
        
#stop(here"); sav_txt_lst <- glb_txt_chr_lst    
    for (txt_var in glbFeatsText) {
        re_str <- "\\b(N|S|E|W|C)( |\\.)(\\w)+"
        print(sprintf("Remaining #%s# terms in %s: ", re_str, txt_var))        
        txt_vctr <- glb_txt_chr_lst[[txt_var]]                
        print(orderBy(~ -.n +pattern, subset(chk_pattern_freq(re_str, ignore.case=FALSE), 
                                             grepl(".", pattern))))
        #grep("N Weaver", txt_vctr, value=TRUE, ignore.case=FALSE)        
    }    

    for (txt_var in glbFeatsText) {
        re_str <- "\\b(North|South|East|West|Central)( |\\.)(\\w)+"
        print(sprintf("Remaining #%s# terms in %s: ", re_str, txt_var))        
        txt_vctr <- glb_txt_chr_lst[[txt_var]]                        
        if (nrow(filtered_df <- subset(chk_pattern_freq(re_str, ignore.case=FALSE), 
                                             grepl(".", pattern))) > 0)
            print(orderBy(~ -.n +pattern, filtered_df))
        #grep("Central (African|Bankers|Cast|Italy|Role|Spring)", txt_vctr, value=TRUE, ignore.case=FALSE)
        #grep("East (Africa|Berlin|London|Poland|Rivals|Spring)", txt_vctr, value=TRUE, ignore.case=FALSE)
        #grep("North (American|Korean|West)", txt_vctr, value=TRUE, ignore.case=FALSE)        
        #grep("South (Pacific|Street)", txt_vctr, value=TRUE, ignore.case=FALSE)
        #grep("St\\. Martins", txt_vctr, value=TRUE, ignore.case=FALSE)
    }    

    find_cmpnd_wrds <- function(txt_vctr) {
        # Enhancements:
        #   - arg should be txt_corpus instead of txt_vctr
        
        txt_corpus <- Corpus(VectorSource(txt_vctr))
        txt_corpus <- tm_map(txt_corpus, content_transformer(tolower), lazy = TRUE)
        txt_corpus <- tm_map(txt_corpus, PlainTextDocument, lazy = TRUE)
        txt_corpus <- tm_map(txt_corpus, removePunctuation, 
                             preserve_intra_word_dashes = TRUE, lazy = FALSE)
        
        # Defaulting to Tf since TfIdf with normalize = TRUE throws a warning for empty docs
        terms_mtrx <- as.matrix(TermDocumentMatrix(txt_corpus,
                                                   control = list(weighting = weightTf)))
        terms_df <- orderBy(~ -Tf, data.frame(term = dimnames(terms_mtrx)$Terms,
                                              Tf = rowSums(terms_mtrx)))
        
        cmpnd_df <- subset(terms_df, grepl("-", term))
        if (nrow(cmpnd_df) == 0) {
            print("   No compounded terms found")
            return(FALSE)
        }
        
        txt_compound_filename <- paste0(glb_txt_munge_filenames_pfx, "compound.csv")
        if (!file.exists(txt_compound_filename))
            stop(txt_compound_filename, " not found!")
        filter_df <- read.csv(txt_compound_filename, comment.char="#", strip.white=TRUE)
        cmpnd_df$filter <- FALSE
        for (row_ix in 1:nrow(filter_df))
            cmpnd_df[!cmpnd_df$filter, "filter"] <- 
            grepl(filter_df[row_ix, "rex_str"], 
                  cmpnd_df[!cmpnd_df$filter, "term"], ignore.case=TRUE)
        cmpnd_df <- subset(cmpnd_df, !filter)
        # Bug in tm_map(txt_corpus, removePunctuation, preserve_intra_word_dashes=TRUE) ???
        #   "net-a-porter" gets converted to "net-aporter"
        #grep("net-a-porter", txt_vctr, ignore.case=TRUE, value=TRUE)
        #grep("maser-laser", txt_vctr, ignore.case=TRUE, value=TRUE)
        #txt_corpus[[which(grepl("net-a-porter", txt_vctr, ignore.case=TRUE))]]
        #grep("\\b(across|longer)-(\\w)", cmpnd_df$term, ignore.case=TRUE, value=TRUE)
        #grep("(\\w)-(affected|term)\\b", cmpnd_df$term, ignore.case=TRUE, value=TRUE)
        
        print(sprintf("nrow(cmpnd_df): %d", nrow(cmpnd_df)))
        myprint_df(cmpnd_df)
    }

    # This should be run after glb_txt_corpus_lst is created with tolower
    extract.features_chunk_df <- myadd_chunk(extract.features_chunk_df, 
            paste0("extract.features_", "process.text_reporting_compound_terms"), major.inc=FALSE)
    
#     for (txt_var in glbFeatsText) {
#         print(sprintf("Remaining compound terms in %s: ", txt_var))        
#         find_cmpnd_wrds(txt_vctr = glb_txt_chr_lst[[txt_var]])
#         #grep("thirty-five", txt_vctr, ignore.case=TRUE, value=TRUE)
#         #rex_str <- glb_txt_map_df[grepl("hirty", glb_txt_map_df$rex_str), "rex_str"]
#     }

    extract.features_chunk_df <- myadd_chunk(extract.features_chunk_df, 
            paste0("extract.features_", "build.corpus"), major.inc=TRUE)
    
    get_txt_terms <- function(terms_TDM) {
        terms_mtrx <- as.matrix(as.TermDocumentMatrix(terms_TDM))
        docms_mtrx <- as.matrix(as.DocumentTermMatrix(terms_TDM))        
        terms_df <- data.frame(term = dimnames(terms_mtrx)$Terms,
                               weight = rowSums(terms_mtrx),
                               freq = rowSums(terms_mtrx > 0))
        terms_df$pos <- 1:nrow(terms_df)
        terms_df$cor.y <- 
            cor(docms_mtrx[glbObsAll$.src == "Train",], 
                as.numeric(glbObsAll[glbObsAll$.src == "Train", glb_rsp_var]),
                              use = "pairwise.complete.obs")
        terms_df$cor.y.abs <- abs(terms_df$cor.y)
#         .rnorm.cor.y.abs <- abs(cor(glbObsAll[glbObsAll$.src == "Train", ".rnorm"],
#                         as.numeric(glbObsAll[glbObsAll$.src == "Train", glb_rsp_var]),
#                                 use = "pairwise.complete.obs"))
        
        terms_df$chisq.stat <- NA; terms_df$chisq.pval <- NA
        for (ix in 1:nrow(terms_df)) {
        #for (ix in 1:743) {
            if (length(unique(docms_mtrx[glbObsAll$.src == "Train", ix])) > 1) {
                chisq <- chisq.test(docms_mtrx[glbObsAll$.src == "Train", ix], 
                                    glbObsAll[glbObsAll$.src == "Train", glb_rsp_var])
                terms_df[ix, "chisq.stat"] <- chisq$statistic
                terms_df[ix, "chisq.pval"] <- chisq$p.value 
            }
        }
        
        nzv_df <- nzv(docms_mtrx[glbObsAll$.src == "Train",], freqCut = glb_nzv_freqCut,
                   uniqueCut = glb_nzv_uniqueCut, saveMetrics = TRUE)
        terms_df$nzv.freqRatio <- nzv_df$freqRatio
#         terms_df$nzv.freqRatio.cut.fctr <- cut(terms_df$nzv.freqRatio, 
#                                                 breaks = sort(c(min(terms_df$nzv.freqRatio), 
#                                                                 glb_nzv_freqCut,
#                                                                 max(terms_df$nzv.freqRatio))))
        terms_df$nzv.percentUnique <- nzv_df$percentUnique
#         terms_df$nzv.percentUnique.cut.fctr <- cut(terms_df$nzv.percentUnique, 
#                     breaks = sort(c(min(terms_df$nzv.percentUnique) - .Machine$double.neg.eps, 
#                                                             glb_nzv_uniqueCut,
#                                                             max(terms_df$nzv.percentUnique))))
#         terms_df$nzv.quad.fctr <- as.factor(paste0("fRatio:", terms_df$nzv.freqRatio.cut.fctr,
#                                             "\n%Unq:", terms_df$nzv.percentUnique.cut.fctr))
        terms_df$nzv <- nzv_df$nzv
        
        for (cls in unique(glbObsAll[, glb_txt_cor_var])) {
            if (!is.na(cls))
                terms_df[, paste0("weight.", as.character(cls))] <- 
                    colSums(t(terms_mtrx) * 
                            as.numeric(!is.na(glbObsAll[, glb_txt_cor_var]) &
                                        (glbObsAll[, glb_txt_cor_var] == cls))) else
                terms_df[, paste0("weight.", as.character(cls))] <- 
                    colSums(t(terms_mtrx) * 
                            as.numeric(is.na(glbObsAll[, glb_txt_cor_var])))
        }    
        
        # Check all calls to get_terms_DTM_terms to change returned order assumption
        return(terms_df <- orderBy(~ -weight, terms_df))
    }
    #plt_full_df <- get_terms_DTM_terms(terms_DTM=glb_full_terms_DTM_lst[[txt_var]])
    
    get_corpus_terms <- function(txt_corpus) {
        return(terms_df <- get_txt_terms(terms_TDM = 
                        TermDocumentMatrix(txt_corpus, control = glb_txt_terms_control)))
    }
    
    myreplacePunctuation <- function(x) {
        return(gsub("[[:punct:]]+", " ", gsub("('|-)", "", x)))
    }
    
#stop(here"); glb_to_sav()    
    glb_txt_corpus_lst <- list()
    print(sprintf("Building glb_txt_corpus_lst..."))
    glb_txt_corpus_lst <- foreach(txt_var = glbFeatsText, .verbose = FALSE) %dopar% {
    #glb_txt_corpus_lst <- foreach(txt_var = glbFeatsText, .verbose = TRUE) %do% {        
    #for (txt_var in glbFeatsText) {
        txt_corpus <- Corpus(VectorSource(glb_txt_chr_lst[[txt_var]]))
        txt_corpus <- tm_map(txt_corpus, PlainTextDocument, lazy = TRUE)
        txt_corpus <- tm_map(txt_corpus, content_transformer(tolower), lazy = TRUE) #nuppr
        # removePunctuation does not replace with whitespace. Use a custom transformer ???
        #txt_corpus <- tm_map(txt_corpus, removePunctuation, lazy = TRUE) #npnct<chr_ix>
        txt_corpus <- tm_map(txt_corpus, content_transformer(myreplacePunctuation)
                             , lazy = TRUE) #npnct<chr_ix>
#         txt-corpus <- tm_map(txt_corpus, content_transformer(function(x, pattern) gsub(pattern, "", x))   
        if (!is.null(glb_txt_stop_words[[txt_var]]))
            txt_corpus <- tm_map(txt_corpus, removeWords, glb_txt_stop_words[[txt_var]],
                                 lazy = FALSE)#, lazy=TRUE) #nstopwrds

        # foreach result is based on .Last.Eval
        txt_corpus <- txt_corpus
        # glb_txt_corpus_lst[[txt_var]] <- txt_corpus
    }
    names(glb_txt_corpus_lst) <- glbFeatsText
    
mycombineSynonyms <- content_transformer(function(x, syn=NULL) { 
    Reduce(function(a,b) {
        gsub(paste0("\\b(", paste(b$syns, collapse = "|"),")\\b"), b$word, a)}, syn, x)   
})    
    
#stop(here"); glb_to_sav(); sav_txt_corpus <- glb_txt_corpus_lst[[txt_var]]; all.equal(sav_txt_corpus, glb_txt_corpus_lst[[txt_var]]); glb_txt_corpus_lst[[txt_var]] <- sav_txt_corpus
    glb_post_stop_words_terms_df_lst <- list(); 
    glb_post_stop_words_terms_mtrx_lst <- list();     
    glb_post_stem_words_terms_df_lst <- list(); 
    glb_post_stem_words_terms_mtrx_lst <- list();     
    for (txt_var in glbFeatsText) {
        print(sprintf("    Top_n post-stop term weights for %s:", txt_var))
        # This impacts stemming probably due to lazy parameter
        print(myprint_df(full_terms_df <-
                             get_corpus_terms(txt_corpus=glb_txt_corpus_lst[[txt_var]]), 
                        glbFeatsTextTermsMax[[txt_var]]))
        glb_post_stop_words_terms_df_lst[[txt_var]] <- full_terms_df
        terms_stop_mtrx <- as.matrix(DocumentTermMatrix(glb_txt_corpus_lst[[txt_var]], 
                                        control=glb_txt_terms_control))
        rownames(terms_stop_mtrx) <- rownames(glbObsAll) # print undreadable otherwise
        glb_post_stop_words_terms_mtrx_lst[[txt_var]] <- terms_stop_mtrx
        
        tmp_allobs_df <- glbObsAll[, c(glb_id_var, glb_rsp_var)]
        tmp_allobs_df$terms.post.stop.n <- rowSums(terms_stop_mtrx > 0)
        tmp_allobs_df$terms.post.stop.n.log <- log(1 + tmp_allobs_df$terms.post.stop.n)
        tmp_allobs_df$weight.post.stop.sum <- rowSums(terms_stop_mtrx)        
        
        print(sprintf("    Top_n stem term weights for %s:", txt_var))        
        glb_txt_corpus_lst[[txt_var]] <- tm_map(glb_txt_corpus_lst[[txt_var]], stemDocument,
                                            "english", lazy=FALSE)
        if (!is.null(glb_txt_synonyms[[txt_var]])) {
            syn_lst <- myrmNullObj(glb_txt_synonyms[[txt_var]])
            glb_txt_corpus_lst[[txt_var]] <- tm_map(glb_txt_corpus_lst[[txt_var]],
                                                    mycombineSynonyms,
                                                    syn_lst, lazy=FALSE)
        }    
        
        print(myprint_df(full_terms_df <- get_corpus_terms(glb_txt_corpus_lst[[txt_var]]), 
                   glbFeatsTextTermsMax[[txt_var]]))
        glb_post_stem_words_terms_df_lst[[txt_var]] <- full_terms_df        
        terms_stem_mtrx <- as.matrix(DocumentTermMatrix(glb_txt_corpus_lst[[txt_var]], 
                                        control=glb_txt_terms_control))
        rownames(terms_stem_mtrx) <- rownames(glbObsAll) # print undreadable otherwise
        glb_post_stem_words_terms_mtrx_lst[[txt_var]] <- terms_stem_mtrx
        
        tmp_allobs_df$terms.post.stem.n <- rowSums(terms_stem_mtrx > 0)
        tmp_allobs_df$terms.post.stem.n.log <- log(1 + tmp_allobs_df$terms.post.stem.n)
        tmp_allobs_df$weight.post.stem.sum <- rowSums(terms_stem_mtrx)
        
        tmp_allobs_df$terms.n.stem.stop.Ratio <- 
            1.0 * tmp_allobs_df$terms.post.stem.n / tmp_allobs_df$terms.post.stop.n
        tmp_allobs_df[(is.nan(tmp_allobs_df$terms.n.stem.stop.Ratio) | 
                       is.infinite(tmp_allobs_df$terms.n.stem.stop.Ratio)), 
                      "terms.n.stem.stop.Ratio"] <- 1.0                
        if ((n.errors <- sum(tmp_allobs_df$terms.n.stem.stop.Ratio > 1)) > 0)
            stop(n.errors, " obs in tmp_allobs_df have terms.n.stem.stop.Ratio > 1", 
                 " happening due to terms filtered by glb_txt_terms_control$bounds$global[1] but stemmable to other terms")
        #print(head(subset(tmp_allobs_df, terms.n.stem.stop.Ratio > 1)))
        #glbObsAll[(row_ix <- which(glbObsAll$UniqueID == 10465)), ]
        #terms_stop_mtrx[row_ix, terms_stop_mtrx[row_ix, ] > 0]
        #setdiff(names(terms_stem_mtrx[row_ix, terms_stem_mtrx[row_ix, ] > 0]), names(terms_stop_mtrx[row_ix, terms_stop_mtrx[row_ix, ] > 0]))
        #mydspObs(list(descr.my.contains="updat"))
        
        tmp_allobs_df$weight.sum.stem.stop.Ratio <- 
            1.0 * tmp_allobs_df$weight.post.stem.sum / tmp_allobs_df$weight.post.stop.sum
        tmp_allobs_df[is.nan(tmp_allobs_df$weight.sum.stem.stop.Ratio) | 
                      is.infinite(tmp_allobs_df$weight.sum.stem.stop.Ratio), 
                      "weight.sum.stem.stop.Ratio"] <- 1.0                
        
        tmp_trnobs_df <- tmp_allobs_df[!is.na(tmp_allobs_df[, glb_rsp_var]), ]
        print(cor(as.matrix(tmp_trnobs_df[, -c(1, 2)]), 
                  as.numeric(tmp_trnobs_df[, glb_rsp_var])))

        txt_var_pfx <- toupper(substr(txt_var, 1, 1))
        tmp_allobs_df <- tmp_allobs_df[, -c(1, 2)]
        names(tmp_allobs_df) <- paste(paste0(txt_var_pfx, "."), names(tmp_allobs_df), sep = "")
        glbObsAll <- cbind(glbObsAll, tmp_allobs_df)
        glbFeatsExclude <- c(glbFeatsExclude, 
                paste(paste0(txt_var_pfx, ".terms.post."), c("stop.n", "stem.n"), sep = ""))
    }
    
    require(wordcloud)
    for (txt_var in glbFeatsText) {
        print(sprintf("    Wordcloud post-stem term weights for %s:", txt_var))
        m <- glb_post_stem_words_terms_mtrx_lst[[txt_var]]
        # calculate the frequency of words
        v <- sort(colSums(m), decreasing = TRUE)
        myNames <- names(v)
        d <- data.frame(word = myNames, freq = v)
        wordcloud(d$word, d$freq, min.freq = glb_txt_terms_control$bounds$global[1])
    }    

    for (txt_var in glbFeatsText) {    
        .rnorm.cor.y.abs <- abs(cor(glbObsAll[glbObsAll$.src == "Train", ".rnorm"],
                    as.numeric(glbObsAll[glbObsAll$.src == "Train", glb_rsp_var]),
                                use = "pairwise.complete.obs"))
        plt_df <- subset(glb_post_stem_words_terms_df_lst[[txt_var]], !is.na(cor.y))
        plt_df$nzv.freqRatio.cut.fctr <- cut(plt_df$nzv.freqRatio, 
                                            breaks = sort(c(min(plt_df$nzv.freqRatio), 
                                                                glb_nzv_freqCut,
                                                            max(plt_df$nzv.freqRatio))))
        plt_df$nzv.percentUnique.cut.fctr <- cut(plt_df$nzv.percentUnique, 
                breaks = sort(c(min(plt_df$nzv.percentUnique) - .Machine$double.neg.eps, 
                                                            glb_nzv_uniqueCut,
                                                        max(plt_df$nzv.percentUnique))))
        plt_df$nzv.quad.fctr <- as.factor(paste0("fRatio:", plt_df$nzv.freqRatio.cut.fctr,
                                            "\n%Unq:", plt_df$nzv.percentUnique.cut.fctr))
        labelCnd <- !plt_df$nzv | 
                    (!is.na(plt_df$chisq.pval) & (plt_df$chisq.pval < 0.05)) & 
                (!is.na(plt_df$cor.y.abs) & (plt_df$cor.y.abs > .rnorm.cor.y.abs))
        plt_df$label <- NA; plt_df[labelCnd, "label"] <- plt_df[labelCnd, "term"]
        print(ggplot(plt_df, aes(x = cor.y, y = chisq.stat)) + 
                  geom_point(aes(color = nzv.quad.fctr, size = weight)) +
                  geom_text(aes(label = label), color = "gray50") + 
                  # geom_vline(xintercept = 0) + 
            geom_vline(xintercept = c(-1, +1) * .rnorm.cor.y.abs, color = "gray", 
                       linetype = "dashed") + 
            ggtitle(txt_var))
    }
        
    extract.features_chunk_df <- myadd_chunk(extract.features_chunk_df, 
            paste0("extract.features_", "extract.DTM"), major.inc=TRUE)

    glb_full_DTM_lst <- list(); glb_sprs_DTM_lst <- list();
    for (txt_var in glbFeatsText) {
        print(sprintf("Extracting term weights for %s...", txt_var))        
        txt_corpus <- glb_txt_corpus_lst[[txt_var]]
        
        full_DTM <- DocumentTermMatrix(txt_corpus, 
                                          control=glb_txt_terms_control)
        sprs_DTM <- removeSparseTerms(full_DTM, 
                                            glb_sprs_thresholds[txt_var])
        
        glb_full_DTM_lst[[txt_var]] <- full_DTM
        glb_sprs_DTM_lst[[txt_var]] <- sprs_DTM
    }

    extract.features_chunk_df <- myadd_chunk(extract.features_chunk_df, 
            paste0("extract.features_", "report.DTM"), major.inc=TRUE)
    
    require(reshape2)
    for (txt_var in glbFeatsText) {
        print(sprintf("Reporting term weights for %s...", txt_var))        
        full_DTM <- glb_full_DTM_lst[[txt_var]]
        sprs_DTM <- glb_sprs_DTM_lst[[txt_var]]        

        print("   Full TermMatrix:"); print(full_DTM)
        full_terms_df <- get_txt_terms(full_DTM)
#         full_terms_df <- full_terms_df[, c(2, 1, 3, 4)]
#         col_names <- names(full_terms_df)
#         col_names[2:length(col_names)] <- 
#             paste(col_names[2:length(col_names)], ".full", sep="")
#         names(full_terms_df) <- col_names

        print("   Sparse TermMatrix:"); print(sprs_DTM)
        sprs_terms_df <- get_txt_terms(sprs_DTM)
#         sprs_terms_df <- sprs_terms_df[, c(2, 1, 3, 4)]
#         col_names <- names(sprs_terms_df)
#         col_names[2:length(col_names)] <- 
#             paste(col_names[2:length(col_names)], ".sprs", sep="")
#         names(sprs_terms_df) <- col_names

        #intersect(names(full_terms_df), names(sprs_terms_df))
        terms_df <- merge(full_terms_df, sprs_terms_df, by = c("term", "weight", "freq",
                                        grep("weight\\.", names(full_terms_df), value = TRUE)),
                          all.x = TRUE, suffixes = c(".full", ".sprs"))
        terms_df$in.sprs <- !is.na(terms_df$pos.sprs)
        plt_terms_df <- subset(terms_df, 
                        weight >= min(terms_df$weight[!is.na(terms_df$pos.sprs)], na.rm=TRUE))
        plt_terms_df$label <- ""
        plt_terms_df[is.na(plt_terms_df$pos.sprs), "label"] <- 
            plt_terms_df[is.na(plt_terms_df$pos.sprs), "term"]
#         glb_important_terms[[txt_var]] <- union(glb_important_terms[[txt_var]],
#             plt_terms_df[is.na(plt_terms_df$TfIdf.sprs), "term"])
        print(myplot_scatter(plt_terms_df, "freq", "weight", 
                             colorcol_name="in.sprs") + 
                  geom_text(aes(label=label), color="Black", size=3.5))
        
        melt_terms_df <- orderBy(~ -value, 
                            melt(terms_df, id.vars="term", measure.vars = c("weight", "freq")))
        print(ggplot(melt_terms_df, aes(value, color=variable)) + stat_ecdf() + 
                  geom_hline(yintercept=glb_sprs_thresholds[txt_var], 
                             linetype = "dotted"))
        
        melt_terms_df <- orderBy(~ -value, 
                        melt(subset(terms_df, in.sprs), id.vars="term",
                             measure.vars=grep("weight.", names(terms_df), value=TRUE)))
        print(myplot_hbar(melt_terms_df, "term", "value", colorcol_name="variable"))
        
        melt_terms_df <- orderBy(~ -value, 
                        melt(subset(terms_df, !in.sprs), id.vars="term",
                             measure.vars=grep("weight.", names(terms_df), value=TRUE)))
        print(myplot_hbar(head(melt_terms_df, glbFeatsTextTermsMax[[txt_var]]), "term", "value",
                          colorcol_name="variable"))
    }

#     sav_full_DTM_lst <- glb_full_DTM_lst
#     print(identical(sav_glb_txt_corpus_lst, glb_txt_corpus_lst))
#     print(all.equal(length(sav_glb_txt_corpus_lst), length(glb_txt_corpus_lst)))
#     print(all.equal(names(sav_glb_txt_corpus_lst), names(glb_txt_corpus_lst)))
#     print(all.equal(sav_glb_txt_corpus_lst[["Headline"]], glb_txt_corpus_lst[["Headline"]]))

#     print(identical(sav_full_DTM_lst, glb_full_DTM_lst))
        
    rm(full_terms_mtrx)

    # Create txt features
    if ((length(glbFeatsText) > 1) &&
        (length(unique(pfxs <- sapply(glbFeatsText, 
                    function(txt) toupper(substr(txt, 1, 1))))) < length(glbFeatsText)))
            stop("Prefixes for corpus freq terms not unique: ", pfxs)
    
    extract.features_chunk_df <- myadd_chunk(extract.features_chunk_df, 
                            paste0("extract.features_", "bind.DTM"), 
                                         major.inc=TRUE)
#stop(here"); glb_to_sav(); all.equal(sav_allobs_df, glbObsAll); glbObsAll <- sav_allobs_df
    require(tidyr)
    for (txt_var in glbFeatsText) {
        print(sprintf("Binding DTM for %s...", txt_var))
        txt_var_pfx <- toupper(substr(txt_var, 1, 1))
        
        txt_full_X_df <- as.data.frame(as.matrix(glb_full_DTM_lst[[txt_var]]))
        terms_full_df <- get_txt_terms(glb_full_DTM_lst[[txt_var]])     
        # make.names adds a period to R keywords e.g. "in", "function"
        colnames(txt_full_X_df) <- paste(txt_var_pfx, ".T.",
                                    make.names(colnames(txt_full_X_df)), sep="")
        rownames(txt_full_X_df) <- rownames(glbObsAll) # warning otherwise
        
#         plt_full_df <- terms_full_df
#         names(plt_full_df)[grepl("weight$", names(plt_full_df))] <- "weight.all"
#     #     gather(plt_full_df[1:5, ], domain, TfIdf, -matches("!(TfIdf)"))
#     #     gather(plt_full_df[1:5, grepl("TfIdf", names(plt_full_df))], domain, TfIdf) 
#     #     gather(plt_full_df[1:5, ], domain, TfIdf, 
#     #            -names(plt_full_df)[!grepl("TfIdf", names(plt_full_df))]) 
#         plt_full_df <- gather(plt_full_df, domain, weight, 
#                               -c(term, freq, pos, cor.y, cor.y.abs))
#         plt_full_df$label <- NA
#         top_val_terms <- orderBy(~-weight, terms_full_df)$term[1:glbFeatsTextTermsMax[[txt_var]]]
#         plt_full_df[plt_full_df$term %in% top_val_terms, "label"] <- 
#             plt_full_df[plt_full_df$term %in% top_val_terms, "term"]
#         top_cor_terms <- orderBy(~-cor.y.abs,
#                                  terms_full_df)$term[1:glbFeatsTextTermsMax[[txt_var]]]
#         plt_full_df[plt_full_df$term %in% top_cor_terms, "label"] <- 
#             plt_full_df[plt_full_df$term %in% top_cor_terms, "term"]
#         #plt_full_df$type <- "none"
#         plt_full_df[plt_full_df$term %in% top_val_terms, "type"] <- "top.weight" 
#         plt_full_df[plt_full_df$term %in% top_cor_terms, "type"] <- "top.cor"
#         plt_full_df[plt_full_df$term %in% intersect(top_val_terms, top_cor_terms), "type"] <-
#             "top.both"
#         cor.y.rnorm <- cor(glbObsAll$.rnorm, as.numeric(glbObsAll[, glb_rsp_var]),
#                            use = "pairwise.complete.obs")
#         print(ggplot(plt_full_df, aes(x=weight, y=cor.y)) + facet_wrap(~ domain) + 
#                 geom_point(aes(size=freq), color="grey") + 
#                 geom_jitter() + 
#                 geom_text(aes(label=label, color=type), size=3.5) +
#         #geom_hline(yintercept=cor.y.rnorm, color="red") + 
#         geom_hline(yintercept=c(cor.y.rnorm, -cor.y.rnorm), color="red"))

#stop(here"); glb_to_sav()        
        if (glbFeatsTextFilter == "sparse") {
            txt_X_df <- as.data.frame(as.matrix(glb_sprs_DTM_lst[[txt_var]]))
            select_terms <- make.names(colnames(txt_X_df))
#             colnames(txt_X_df) <- paste(txt_var_pfx, ".T.",
#                                         make.names(colnames(txt_X_df)), sep="")
#             rownames(txt_X_df) <- rownames(glbObsAll) # warning otherwise
        } else if (glbFeatsTextFilter == "top.val") {
            select_terms <- orderBy(~-weight,
                                    terms_full_df)$term[1:glbFeatsTextTermsMax[[txt_var]]]
#             txt_X_df <- txt_full_X_df[, subset(terms_full_df, term %in% select_terms)$pos,
#                                       FALSE]
        } else if (glbFeatsTextFilter == "top.cor") {
            select_terms <- orderBy(~-cor.y.abs,
                                    terms_full_df)$term[1:glbFeatsTextTermsMax[[txt_var]]]
#             txt_X_df <- txt_full_X_df[, subset(terms_full_df, term %in% select_terms)$pos,
#                                       FALSE]
        } else if (glbFeatsTextFilter == "top.chisq") {
            select_terms <- orderBy(~-chisq.stat,
                                    subset(terms_full_df, chisq.pval < 0.05)
                                    )$term[1:glbFeatsTextTermsMax[[txt_var]]]
        } else if (glbFeatsTextFilter == "union.top.val.cor") {
            select_terms <- union(
                orderBy(~-weight   , terms_full_df)$term[1:glbFeatsTextTermsMax[[txt_var]]],
                orderBy(~-cor.y.abs, terms_full_df)$term[1:glbFeatsTextTermsMax[[txt_var]]])
        } else stop(
        "glbFeatsTextFilter should be one of c('sparse', 'top.val', 'top.cor', 'union.top.val.cor', 'top.chisq') vs. '",
                    glbFeatsTextFilter, "'")    
        
        assoc_terms_lst <- findAssocs(glb_full_DTM_lst[[txt_var]], select_terms, 
                                      glbFeatsTextAssocCor[[txt_var]])
        assoc_terms <- c(NULL)
        for (term in names(assoc_terms_lst))
            if (length(assoc_terms_lst[[term]]) > 0)
                assoc_terms <- union(assoc_terms, names(assoc_terms_lst[[term]]))

#stop(here"); glb_to_sav()
        txt_X_df <- txt_full_X_df[, 
                        subset(terms_full_df, term %in% c(select_terms, assoc_terms))$pos,
                                    FALSE]
        glbObsAll <- cbind(glbObsAll, txt_X_df) # TfIdf is normalized
        #glbObsAll <- cbind(glbObsAll, log_X_df) # if using non-normalized metrics 
    }
    #identical(chk_entity_df, glbObsAll)
    #chk_entity_df <- glbObsAll

    extract.features_chunk_df <- myadd_chunk(extract.features_chunk_df, 
                            paste0("extract.features_", "bind.DXM"), 
                                         major.inc=TRUE)

#stop(here"); sav_allobs_df <- glbObsAll; glbObsAll <- sav_allobs_df
    glb_punct_vctr <- c("!", "\"", "#", "\\$", "%", "&", "'", 
                        "\\(|\\)",# "\\(", "\\)", 
                        "\\*", "\\+", ",", "-", "\\.", "/", ":", ";", 
                        "<|>", # "<", 
                        "=", 
                        # ">", 
                        "\\?", "@", "\\[", "\\\\", "\\]", "\\^", "_", "`", 
                        "\\{", "\\|", "\\}", "~")
    txt_X_df <- glbObsAll[, c(glb_id_var, ".rnorm"), FALSE]
    txt_X_df <- foreach(txt_var=glbFeatsText, .combine=cbind) %dopar% {   
    #for (txt_var in glbFeatsText) {
        print(sprintf("Binding DXM for %s...", txt_var))
        txt_var_pfx <- toupper(substr(txt_var, 1, 1))        
        
        txt_full_DTM_mtrx <- as.matrix(glb_full_DTM_lst[[txt_var]])
        rownames(txt_full_DTM_mtrx) <- rownames(glbObsAll) # print undreadable otherwise
        #print(txt_full_DTM_mtrx[txt_full_DTM_mtrx[, "ebola"] != 0, "ebola"])
        
        # Create <txt_var>.T.<term> for glb_important_terms
        for (term in glb_important_terms[[txt_var]])
            txt_X_df[, paste0(txt_var_pfx, ".T.", make.names(term))] <- 
                txt_full_DTM_mtrx[, term]
                
        # Create <txt_var>.wrds.n.log & .wrds.unq.n.log
        txt_X_df[, paste0(txt_var_pfx, ".wrds.n.log")] <- 
            log(1 + mycount_pattern_occ("\\w+", glb_txt_chr_lst[[txt_var]]))
        txt_X_df[, paste0(txt_var_pfx, ".wrds.unq.n.log")] <- 
            log(1 + rowSums(txt_full_DTM_mtrx != 0))
        txt_X_df[, paste0(txt_var_pfx, ".weight.sum")] <- 
            rowSums(txt_full_DTM_mtrx) 
        txt_X_df[, paste0(txt_var_pfx, ".ratio.weight.sum.wrds.n")] <- 
            txt_X_df[, paste0(txt_var_pfx, ".weight.sum")] / 
            (exp(txt_X_df[, paste0(txt_var_pfx, ".wrds.n.log")]) - 1)
        txt_X_df[is.nan(txt_X_df[, paste0(txt_var_pfx, ".ratio.weight.sum.wrds.n")]),
                 paste0(txt_var_pfx, ".ratio.weight.sum.wrds.n")] <- 0

        # Create <txt_var>.chrs.n.log
        txt_X_df[, paste0(txt_var_pfx, ".chrs.n.log")] <- 
            log(1 + mycount_pattern_occ(".", glbObsAll[, txt_var]))
        txt_X_df[, paste0(txt_var_pfx, ".chrs.uppr.n.log")] <- 
            log(1 + mycount_pattern_occ("[[:upper:]]", glbObsAll[, txt_var]))
        txt_X_df[, paste0(txt_var_pfx, ".dgts.n.log")] <- 
            log(1 + mycount_pattern_occ("[[:digit:]]", glbObsAll[, txt_var]))

        # Create <txt_var>.npnct?.log
        # would this be faster if it's iterated over each row instead of 
        #   each created column ???
        for (punct_ix in 1:length(glb_punct_vctr)) { 
#             smp0 <- " "
#             smp1 <- "! \" # $ % & ' ( ) * + , - . / : ; < = > ? @ [ \ ] ^ _ ` { | } ~"
#             smp2 <- paste(smp1, smp1, sep=" ")
#             print(sprintf("Testing %s pattern:", glb_punct_vctr[punct_ix])) 
#             results <- mycount_pattern_occ(glb_punct_vctr[punct_ix], c(smp0, smp1, smp2))
#             names(results) <- NULL; print(results)
            txt_X_df[, 
                paste0(txt_var_pfx, ".chrs.pnct", sprintf("%02d", punct_ix), ".n.log")] <-
                log(1 + mycount_pattern_occ(glb_punct_vctr[punct_ix], 
                                            glbObsAll[, txt_var]))
        }
#         print(head(glbObsAll[glbObsAll[, "A.npnct23.log"] > 0, 
#                                     c("UniqueID", "Popular", "Abstract", "A.npnct23.log")]))    
        
        # Create <txt_var>.niso8859.log
        txt_X_df[, paste0(txt_var_pfx, ".chrs.iso8859.n.log")] <- 
            log1p(mycount_pattern_occ("&#[[:digit:]]{3};", glbObsAll[, txt_var]))
                
        # Create <txt_var>.wrds.stop.n.log & <txt_var>ratio.wrds.stop.n.wrds.n
        if (!is.null(glb_txt_stop_words[[txt_var]])) {
            stop_words_rex_str <- paste0("\\b(", 
                                         paste0(glb_txt_stop_words[[txt_var]], collapse="|"),
                                         ")\\b")
            txt_X_df[, paste0(txt_var_pfx, ".wrds.stop.n", ".log")] <-
                log(1 + mycount_pattern_occ(stop_words_rex_str, glb_txt_chr_lst[[txt_var]]))
            txt_X_df[, paste0(txt_var_pfx, ".ratio.wrds.stop.n.wrds.n")] <-
                exp(txt_X_df[, paste0(txt_var_pfx, ".wrds.stop.n", ".log")] - 
                    txt_X_df[, paste0(txt_var_pfx, ".wrds.n", ".log")])
        }

        # Create <txt_var>.P.http
        txt_X_df[, paste(txt_var_pfx, ".P.http", sep="")] <- 
            log1p(mycount_pattern_occ("http", glbObsAll[, txt_var]))
            
        # Create <txt_var>.P.<user-spec-pattern>
        for (pattern in names(glbFeatsTextPatterns[[txt_var]]))
            txt_X_df[, paste(txt_var_pfx, ".P.", pattern, sep="")] <- 
                log1p(mycount_pattern_occ(glbFeatsTextPatterns[[txt_var]][pattern], 
                                            glbObsAll[, txt_var]))
    
        txt_X_df <- subset(txt_X_df, select=-.rnorm)
        txt_X_df <- txt_X_df[, -grep(glb_id_var, names(txt_X_df), fixed=TRUE), FALSE]
        #glbObsAll <- cbind(glbObsAll, txt_X_df)
    }
    glbObsAll <- cbind(glbObsAll, txt_X_df)
    #myplot_box(glbObsAll, "A.sum.TfIdf", glb_rsp_var)

    # Generate summaries
#     print(summary(glbObsAll))
#     print(sapply(names(glbObsAll), function(col) sum(is.na(glbObsAll[, col]))))
#     print(summary(glbObsTrn))
#     print(sapply(names(glbObsTrn), function(col) sum(is.na(glbObsTrn[, col]))))
#     print(summary(glbObsNew))
#     print(sapply(names(glbObsNew), function(col) sum(is.na(glbObsNew[, col]))))

    glbFeatsExclude <- union(glbFeatsExclude, 
                                          glbFeatsText)
    rm(log_X_df, txt_X_df)
}

# print(sapply(names(glbObsTrn), function(col) sum(is.na(glbObsTrn[, col]))))
# print(sapply(names(glbObsNew), function(col) sum(is.na(glbObsNew[, col]))))

# print(myplot_scatter(glbObsTrn, "<col1_name>", "<col2_name>", smooth=TRUE))

#stop(here"); glb_to_sav(); glbObsAll <- sav_allobs_df
if (!is.null(glbFeatsPrice)) {
    for (var in glbFeatsPrice) {
        for (digit in 1:(log10(max(glbObsAll[, var], na.rm=TRUE)) + 1)) {
            glbObsAll[, paste0(var, ".dgt", digit, ".is9")] <- 
                as.numeric(as.integer((as.integer(glbObsAll[, var]) %% (10 ^ digit)) / 
                                          (10 ^ (digit - 1))) == 9)
#             glbObsAll[, paste0(var, ".dgt", digit, ".is9.fctr")] <- 
#                 as.factor(as.integer((as.integer(glbObsAll[, var]) %% (10 ^ digit)) / 
#                                           (10 ^ (digit - 1))) == 9)
        }
        for (decimal in 1:2) {
            glbObsAll[, paste0(var, ".dcm", decimal, ".is9")] <- 
                as.numeric(as.integer(glbObsAll[, var] * (10 ^ decimal)) %% 10 == 9)
#             glbObsAll[, paste0(var, ".dcm", decimal, ".is9.fctr")] <- 
#                 as.factor(as.integer(glbObsAll[, var] * (10 ^ decimal)) %% 10 == 9)
        }
    }
    #as.numeric((as.integer(startprice) %% 10) == 9)    
}

rm(corpus_lst
   , glb_sprs_DTM_lst #, glb_full_DTM_lst
   , txt_corpus, txt_vctr)
```

```
## Warning in rm(corpus_lst, glb_sprs_DTM_lst, txt_corpus, txt_vctr): object
## 'corpus_lst' not found
```

```
## Warning in rm(corpus_lst, glb_sprs_DTM_lst, txt_corpus, txt_vctr): object
## 'glb_sprs_DTM_lst' not found
```

```
## Warning in rm(corpus_lst, glb_sprs_DTM_lst, txt_corpus, txt_vctr): object
## 'txt_corpus' not found
```

```
## Warning in rm(corpus_lst, glb_sprs_DTM_lst, txt_corpus, txt_vctr): object
## 'txt_vctr' not found
```

```r
extract.features_chunk_df <- myadd_chunk(extract.features_chunk_df, "extract.features_end", 
                                     major.inc=TRUE)
```

```
##                                 label step_major step_minor label_minor
## 2 extract.features_factorize.str.vars          2          0           0
## 3                extract.features_end          3          0           0
##      bgn    end elapsed
## 2 22.552 22.578   0.026
## 3 22.578     NA      NA
```

```r
myplt_chunk(extract.features_chunk_df)
```

```
##                                 label step_major step_minor label_minor
## 2 extract.features_factorize.str.vars          2          0           0
## 1                extract.features_bgn          1          0           0
##      bgn    end elapsed duration
## 2 22.552 22.578   0.026    0.026
## 1 22.539 22.552   0.013    0.013
## [1] "Total Elapsed Time: 22.578 secs"
```

![](NYTBlogs3_base_files/figure-html/extract.features-1.png) 

```r
# if (glb_save_envir)
#     save(glb_feats_df, 
#          glbObsAll, #glbObsTrn, glbObsFit, glbObsOOB, glbObsNew,
#          file=paste0(glb_out_pfx, "extract_features_dsk.RData"))
# load(paste0(glb_out_pfx, "extract_features_dsk.RData"))

replay.petrisim(pn=glb_analytics_pn, 
    replay.trans=(glb_analytics_avl_objs <- c(glb_analytics_avl_objs, 
        "data.training.all","data.new")), flip_coord=TRUE)
```

```
## time	trans	 "bgn " "fit.data.training.all " "predict.data.new " "end " 
## 0.0000 	multiple enabled transitions:  data.training.all data.new model.selected 	firing:  data.training.all 
## 1.0000 	 1 	 2 1 0 0 
## 1.0000 	multiple enabled transitions:  data.training.all data.new model.selected model.final data.training.all.prediction 	firing:  data.new 
## 2.0000 	 2 	 1 1 1 0
```

![](NYTBlogs3_base_files/figure-html/extract.features-2.png) 

```r
#glb_chunks_df <- myadd_chunk(glb_chunks_df, "manage.missing.data", major.inc=TRUE)
```

### Step `3.0: extract features`

```r
glb_chunks_df <- myadd_chunk(glb_chunks_df, "manage.missing.data", major.inc=FALSE)
```

```
##                 label step_major step_minor label_minor    bgn    end
## 5    extract.features          3          0           0 22.532 24.177
## 6 manage.missing.data          3          1           1 24.178     NA
##   elapsed
## 5   1.645
## 6      NA
```

```r
# If mice crashes with error: Error in get(as.character(FUN), mode = "function", envir = envir) : object 'State' of mode 'function' was not found
#   consider excluding 'State' as a feature

# print(sapply(names(glbObsTrn), function(col) sum(is.na(glbObsTrn[, col]))))
# print(sapply(names(glbObsNew), function(col) sum(is.na(glbObsNew[, col]))))
# glbObsTrn <- na.omit(glbObsTrn)
# glbObsNew <- na.omit(glbObsNew)
# df[is.na(df)] <- 0

mycheck_problem_data(glbObsAll, featsExclude = glbFeatsExclude, 
                     fctrMaxUniqVals = glbFctrMaxUniqVals)
```

```
## [1] "numeric data missing in glbObsAll: "
##      Popular Popular.fctr 
##         1870         1870 
## [1] "numeric data w/ 0s in glbObsAll: "
##       WordCount         Popular WordCount.log1p WordCount.root2 
##             109            5439             109             109 
##  WordCount.nexp 
##            2044 
## [1] "numeric data w/ Infs in glbObsAll: "
## named integer(0)
## [1] "numeric data w/ NaNs in glbObsAll: "
## named integer(0)
## [1] "string data missing in glbObsAll: "
##       NewsDesk    SectionName SubsectionName       Headline        Snippet 
##           2408           2899           6176              0             13 
##       Abstract        PubDate    NDSSName.my 
##             17              0              0
```

```r
# glbObsAll <- na.omit(glbObsAll)

# Not refactored into mydsutils.R since glb_*_df might be reassigned
glb_impute_missing_data <- function() {
    
    require(mice)
    set.seed(glb_mice_complete.seed)
    inp_impent_df <- glbObsAll[, setdiff(names(glbObsAll), 
                                union(glbFeatsExclude, glb_rsp_var))]
    print("Summary before imputation: ")
    print(summary(inp_impent_df))
    out_impent_df <- complete(mice(inp_impent_df))
    print(summary(out_impent_df))
    
    ret_vars <- sapply(names(out_impent_df), 
                       function(col) ifelse(!identical(out_impent_df[, col],
                                                       inp_impent_df[, col]), 
                                            col, ""))
    ret_vars <- ret_vars[ret_vars != ""]
    
    # complete(mice()) changes attributes of factors even though values don't change
    for (col in ret_vars) {
        if (inherits(out_impent_df[, col], "factor")) {
            if (identical(as.numeric(out_impent_df[, col]), 
                          as.numeric(inp_impent_df[, col])))
                ret_vars <- setdiff(ret_vars, col)
        }
    }
    return(out_impent_df[, ret_vars])
}

if (glb_impute_na_data && 
    (length(myfind_numerics_missing(glbObsAll)) > 0) &&
    (ncol(nonna_df <- glb_impute_missing_data()) > 0)) {
    for (col in names(nonna_df)) {
        glbObsAll[, paste0(col, ".nonNA")] <- nonna_df[, col]
        glbFeatsExclude <- c(glbFeatsExclude, col)        
    }
}    
    
mycheck_problem_data(glbObsAll, featsExclude = glbFeatsExclude, 
                     fctrMaxUniqVals = glbFctrMaxUniqVals, terminate = TRUE)
```

```
## [1] "numeric data missing in glbObsAll: "
##      Popular Popular.fctr 
##         1870         1870 
## [1] "numeric data w/ 0s in glbObsAll: "
##       WordCount         Popular WordCount.log1p WordCount.root2 
##             109            5439             109             109 
##  WordCount.nexp 
##            2044 
## [1] "numeric data w/ Infs in glbObsAll: "
## named integer(0)
## [1] "numeric data w/ NaNs in glbObsAll: "
## named integer(0)
## [1] "string data missing in glbObsAll: "
##       NewsDesk    SectionName SubsectionName       Headline        Snippet 
##           2408           2899           6176              0             13 
##       Abstract        PubDate    NDSSName.my 
##             17              0              0
```

## Step `3.1: manage missing data`

```r
glb_chunks_df <- myadd_chunk(glb_chunks_df, "cluster.data", major.inc=FALSE)
```

```
##                 label step_major step_minor label_minor    bgn    end
## 6 manage.missing.data          3          1           1 24.178 24.259
## 7        cluster.data          3          2           2 24.259     NA
##   elapsed
## 6   0.081
## 7      NA
```

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
#stop(here"); glb_to_sav(); glbObsAll <- sav_allobs_df
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
                                             by_var=glb_category_var))
    print(sprintf("glbObsAll$%s Entropy: %0.4f (%0.4f pct)",
                    glb_category_var,
            category_ent <- weighted.mean(category_df$.entropy, category_df$.knt),
                    100 * category_ent / allobs_ent))

    glbObsAll$.clusterid <- 1    
    #print(max(table(glbObsAll$myCategory.fctr) / 20))
        
#stop(here"); glb_to_sav()    
    grp_ids <- sort(unique(glbObsAll[, glb_category_var]))
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
        ctgry_allobs_df <- glbObsAll[glbObsAll[, glb_category_var] == grp, ]
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
            c(glb_id_var, glb_cluster_entropy_var, glb_category_var, glbFeatsText, cluster_vars)])
    
        min_dstns_mtrx <- dstns_mtrx
        diag(min_dstns_mtrx) <- 1
        # Float representations issue -2.22e-16 vs. 0.0000
        print(sprintf("min distance(%0.4f) pair:", min(min_dstns_mtrx)))
        row_ix <- ceiling(which.min(min_dstns_mtrx) / ncol(min_dstns_mtrx))
        col_ix <- which.min(min_dstns_mtrx[row_ix, ])
        print(ctgry_allobs_df[c(row_ix, col_ix), 
            c(glb_id_var, glb_cluster_entropy_var, glb_category_var, glbFeatsText,
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
        glbObsAll[glbObsAll[, glb_category_var] == grp,]$.clusterid <-
            clusterGroups
        
        pltIx <- pltIx + 1
    }
    dev.off()
    #all.equal(sav_allobs_df_clusterid, glbObsAll$.clusterid)
    
    print(cluster_df <- mycompute_entropy_df(obs_df=glbObsAll,
                                             entropy_var=glb_cluster_entropy_var,
                                             by_var=glb_category_var))
    print(sprintf("glbObsAll$%s$.clusterid Entropy: %0.4f (%0.4f pct)",
                    glb_category_var,
                cluster_ent <- weighted.mean(cluster_df$.entropy, cluster_df$.knt),
                    100 * cluster_ent / category_ent))
    
    glbObsAll$.clusterid.fctr <- as.factor(glbObsAll$.clusterid)
    # .clusterid.fctr is created automatically (probably ?) later
    glbFeatsExclude <- c(glbFeatsExclude, ".clusterid")
    if (!is.null(glb_category_var))
#         glbFeatsInteractionOnly[ifelse(grepl("\\.fctr", glb_category_var),
#                                             glb_category_var, 
#                                             paste0(glb_category_var, ".fctr"))] <-
#             c(".clusterid.fctr")
        glbFeatsInteractionOnly[[".clusterid.fctr"]] <-
            ifelse(grepl("\\.fctr", glb_category_var), glb_category_var, 
                                                        paste0(glb_category_var, ".fctr"))
            
    if (glbFeatsTextClusterVarsExclude)
        glbFeatsExclude <- c(glbFeatsExclude, cluster_vars)
}

# Last call for data modifications 
#stop(here") # sav_allobs_df <- glbObsAll
# glbObsAll[(glbObsAll$PropR == 0.75) & (glbObsAll$State == "Hawaii"), "PropR.fctr"] <- "N"

# Re-partition
glbObsTrn <- subset(glbObsAll, .src == "Train")
glbObsNew <- subset(glbObsAll, .src == "Test")

glb_chunks_df <- myadd_chunk(glb_chunks_df, "partition.data.training", major.inc=TRUE)
```

```
##                     label step_major step_minor label_minor    bgn    end
## 7            cluster.data          3          2           2 24.259 24.288
## 8 partition.data.training          4          0           0 24.288     NA
##   elapsed
## 7   0.029
## 8      NA
```

## Step `4.0: partition data training`

```r
if (all(is.na(glbObsNew[, glb_rsp_var]))) {
    set.seed(glb_split_sample.seed)
    OOB_size <- nrow(glbObsNew) * 1.1
    
    if (is.null(glb_category_var)) {
        require(caTools)
        split <- sample.split(glbObsTrn[, glb_rsp_var_raw], 
                              SplitRatio=OOB_size / nrow(glbObsTrn))
        glbObsOOB <- glbObsTrn[split ,]            
        glbObsFit <- glbObsTrn[!split, ] 
    } else {
        sample_vars <- c(glb_category_var, glb_rsp_var_raw)
        rspvar_freq_df <- orderBy(reformulate(glb_rsp_var_raw), 
                                 mycreate_sqlxtab_df(glbObsTrn, glb_rsp_var_raw))
        OOB_rspvar_size <- 
            1.0 * OOB_size * rspvar_freq_df$.n / sum(rspvar_freq_df$.n) 
        names(OOB_rspvar_size) <- as.character(rspvar_freq_df[, glb_rsp_var_raw])
        newobs_freq_df <- orderBy(reformulate(glb_category_var),
                                mycreate_sqlxtab_df(glbObsNew, glb_category_var))
        trnobs_freq_df <- orderBy(reformulate(glb_category_var),
                                mycreate_sqlxtab_df(glbObsTrn, glb_category_var))
        ctgry_freq_df <- merge(newobs_freq_df, trnobs_freq_df, 
                                by = glb_category_var,
                                all = TRUE, sort = TRUE, 
                                suffixes = c(".tst", ".trn"))
        ctgry_freq_df[is.na(ctgry_freq_df)] <- 0
        
        obs_freq_df <- mycreate_xtab_df(glbObsTrn, 
                                        c(glb_category_var, glb_rsp_var_raw))
        newobs_freq_df <- orderBy(reformulate(glb_category_var),
                                mycreate_sqlxtab_df(glbObsNew, glb_category_var))
        names(newobs_freq_df) <- gsub(".n", ".n.tst", names(newobs_freq_df), 
                                      fixed = TRUE)
        obs_freq_df <- merge(obs_freq_df, newobs_freq_df, all = TRUE)
        strata_mtrx <- ceiling(
            matrix(obs_freq_df$.n.tst * 1.0 / sum(obs_freq_df$.n.tst)) %*%
                      matrix(OOB_rspvar_size, nrow = 1))
        dimnames(strata_mtrx)[[1]] <- obs_freq_df[, glb_category_var]
        dimnames(strata_mtrx)[[2]] <- 
            as.character(rspvar_freq_df[, glb_rsp_var_raw])
        for (val in rspvar_freq_df[, glb_rsp_var_raw]) {
            trn <- paste0(glb_rsp_var_raw, ".", as.character(val))
            strata <- paste0(".strata.", as.character(val))            
            obs_freq_df[, strata] <- strata_mtrx[, as.character(val)]  
            
            if (length(ix <- which(is.na(obs_freq_df[, trn]))) > 0) {
                # NA obs in a particular category 
                print("Prediction Hints by Catgeory:")
                print(obs_freq_df[ix, ])
                obs_freq_df[ix, strata] <- NA
            }

            if (length((ix <- which(obs_freq_df[, trn] < 
                                    2.0 * obs_freq_df[, strata]))) > 0)
                # More obs in OOB compared to fit currently
                obs_freq_df[ix, strata] <- floor(obs_freq_df[ix, trn] / 2.0)
        
            if (length((ix <- which(obs_freq_df[, strata] == 0))) > 0)
                obs_freq_df[ix, strata] <- 1
        }    
        #print(colSums(obs_freq_df[, -1]))
    
                
        OOB_strata_size <- as.vector(as.matrix(obs_freq_df[, 
                                    grepl("^\\.strata\\.", names(obs_freq_df))]))
        tmp_trnobs_df <- orderBy(reformulate(c(glb_rsp_var_raw, glb_category_var)),
                                glbObsTrn)
        require(sampling)
        split_strata <- sampling::strata(tmp_trnobs_df, 
                               stratanames = c(glb_rsp_var_raw, glb_category_var),
                               size = OOB_strata_size[!is.na(OOB_strata_size)],
                               method = "srswor")
        glbObsOOB <- getdata(tmp_trnobs_df, split_strata)[, names(glbObsTrn)]
        glbObsFit <- glbObsTrn[!glbObsTrn[, glb_id_var] %in% 
                                        glbObsOOB[, glb_id_var], ]
    }
} else {
    print(sprintf("Newdata contains non-NA data for %s; setting OOB to Newdata", 
                  glb_rsp_var))
    glbObsFit <- glbObsTrn; glbObsOOB <- glbObsNew
}
```

```
## [1] "Prediction Hints by Catgeory:"
##    NDSSName.my.fctr Popular.0 Popular.1 .n.tst .strata.0 .strata.1
## 5   #U.S.#Education       325        NA     89        82        17
## 10        Culture##         1        NA     70         1        13
## 12   Foreign#World#       172        NA     47        44         9
## 21          myOther        38        NA      5         5         1
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

```r
if (!is.null(glb_max_fitobs) && (nrow(glbObsFit) > glb_max_fitobs)) {
    warning("glbObsFit restricted to glb_max_fitobs: ", 
            format(glb_max_fitobs, big.mark = ","))
    org_fitobs_df <- glbObsFit
    glbObsFit <- 
        org_fitobs_df[split <- sample.split(org_fitobs_df[, glb_rsp_var_raw], 
                                            SplitRatio = glb_max_fitobs), ]
    org_fitobs_df <- NULL
}

glbObsAll$.lcn <- ""; glbObsTrn$.lcn <- "";
glbObsAll[glbObsAll[, glb_id_var] %in% 
              glbObsFit[, glb_id_var], ".lcn"] <- "Fit"
glbObsTrn[glbObsTrn[, glb_id_var] %in% 
              glbObsFit[, glb_id_var], ".lcn"] <- "Fit"
glbObsAll[glbObsAll[, glb_id_var] %in% 
              glbObsOOB[, glb_id_var], ".lcn"] <- "OOB"
glbObsTrn[glbObsTrn[, glb_id_var] %in% 
              glbObsOOB[, glb_id_var], ".lcn"] <- "OOB"

dsp_class_dstrb <- function(obs_df, location_var, partition_var) {
    xtab_df <- mycreate_xtab_df(obs_df, c(location_var, partition_var))
    rownames(xtab_df) <- xtab_df[, location_var]
    xtab_df <- xtab_df[, -grepl(location_var, names(xtab_df))]
    print(xtab_df)
    print(xtab_df / rowSums(xtab_df, na.rm=TRUE))    
}    

# Ensure proper splits by glb_rsp_var_raw & user-specified feature for OOB vs. new
if (!is.null(glb_category_var)) {
    if (glb_is_classification)
        dsp_class_dstrb(glbObsAll, ".lcn", glb_rsp_var_raw)
    newobs_ctgry_df <- mycreate_sqlxtab_df(subset(glbObsAll, .src == "Test"), 
                                           glb_category_var)
    names(newobs_ctgry_df) <- 
        gsub(".n", ".n.Tst", names(newobs_ctgry_df), fixed = TRUE)
    OOBobs_ctgry_df <- mycreate_sqlxtab_df(subset(glbObsAll, .lcn == "OOB"), 
                                           glb_category_var)
    names(OOBobs_ctgry_df) <- 
        gsub(".n", ".n.OOB", names(OOBobs_ctgry_df), fixed = TRUE)
    fitobs_ctgry_df <- mycreate_sqlxtab_df(subset(glbObsAll, .lcn == "Fit"), 
                                           glb_category_var)
    names(fitobs_ctgry_df) <- 
        gsub(".n", ".n.Fit", names(fitobs_ctgry_df), fixed = TRUE)
    
#     glb_ctgry_df <- merge(newobs_ctgry_df, OOBobs_ctgry_df, by=glb_category_var
#                           , all=TRUE, suffixes=c(".Tst", ".OOB"))
    glb_ctgry_df <- merge(fitobs_ctgry_df, OOBobs_ctgry_df, by = glb_category_var
                          , all = TRUE)
    glb_ctgry_df <- merge(glb_ctgry_df, newobs_ctgry_df, by = glb_category_var
                          , all = TRUE)
    glb_ctgry_df$.freqRatio.Fit <- 
        glb_ctgry_df$.n.Fit / sum(glb_ctgry_df$.n.Fit, na.rm = TRUE)
    glb_ctgry_df$.freqRatio.OOB <- 
        glb_ctgry_df$.n.OOB / sum(glb_ctgry_df$.n.OOB, na.rm = TRUE)
    glb_ctgry_df$.freqRatio.Tst <- 
        glb_ctgry_df$.n.Tst / sum(glb_ctgry_df$.n.Tst, na.rm = TRUE)
    print(orderBy(~-.freqRatio.Tst-.freqRatio.OOB-.freqRatio.Fit, glb_ctgry_df))
}
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
##                      NDSSName.my.fctr .n.Fit .n.OOB .n.Tst .freqRatio.Fit
## 1                                  ##    913    371    342    0.190049958
## 6       Business#BusinessDay#Dealbook    629    323    304    0.130932556
## 11                      Culture#Arts#    490    185    174    0.101998335
## 15                      OpEd#Opinion#    437     89    164    0.090965862
## 9                Business#Technology#    213    126    114    0.044338052
## 19                           TStyle##    623    101    105    0.129683597
## 5                     #U.S.#Education    243     82     89    0.050582848
## 10                          Culture##     NA      1     70             NA
## 14                 Metro#N.Y./Region#    128     70     67    0.026644463
## 18                       Styles#U.S.#    127     50     61    0.026436303
## 16                    Science#Health#    148     48     57    0.030807660
## 13          Foreign#World#AsiaPacific    150     53     56    0.031223980
## 2                        #Multimedia#     92     49     52    0.019150708
## 12                     Foreign#World#    128     44     47    0.026644463
## 8          Business#Crosswords/Games#    105     18     42    0.021856786
## 7  Business#BusinessDay#SmallBusiness    100     40     41    0.020815987
## 20                     Travel#Travel#     83     34     35    0.017277269
## 3              #Opinion#RoomForDebate     42     20     20    0.008742714
## 17                    Styles##Fashion    104     15     15    0.021648626
## 4            #Opinion#ThePublicEditor     16      4     10    0.003330558
## 21                            myOther     33      5      5    0.006869276
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

```r
print("glbObsAll: "); print(dim(glbObsAll))
```

```
## [1] "glbObsAll: "
```

```
## [1] 8402   19
```

```r
print("glbObsTrn: "); print(dim(glbObsTrn))
```

```
## [1] "glbObsTrn: "
```

```
## [1] 6532   19
```

```r
print("glbObsFit: "); print(dim(glbObsFit))
```

```
## [1] "glbObsFit: "
```

```
## [1] 4804   18
```

```r
print("glbObsOOB: "); print(dim(glbObsOOB))
```

```
## [1] "glbObsOOB: "
```

```
## [1] 1728   18
```

```r
print("glbObsNew: "); print(dim(glbObsNew))
```

```
## [1] "glbObsNew: "
```

```
## [1] 1870   18
```

```r
# # Does not handle NULL or length(glb_id_var) > 1

if (glb_save_envir)
    save(glbObsAll, #glbObsTrn, glbObsFit, glbObsOOB, glbObsNew,
         file=paste0(glb_out_pfx, "blddfs_dsk.RData"))
# load(paste0(glb_out_pfx, "blddfs_dsk.RData"))

rm(split)
```

```
## Warning in rm(split): object 'split' not found
```

```r
glb_chunks_df <- myadd_chunk(glb_chunks_df, "select.features", major.inc=TRUE)
```

```
##                     label step_major step_minor label_minor    bgn    end
## 8 partition.data.training          4          0           0 24.288 25.428
## 9         select.features          5          0           0 25.428     NA
##   elapsed
## 8    1.14
## 9      NA
```

## Step `5.0: select features`

```r
#stop(here"); glb_to_sav(); glbObsAll <- sav_allobs_df
print(glb_feats_df <- myselect_features(entity_df=glbObsTrn, 
                       exclude_vars_as_features=glbFeatsExclude, 
                       rsp_var=glb_rsp_var))
```

```
##                                id        cor.y exclude.as.feat   cor.y.abs
## Popular                   Popular  1.000000000               1 1.000000000
## WordCount.root2   WordCount.root2  0.292120679               0 0.292120679
## WordCount               WordCount  0.257526549               1 0.257526549
## WordCount.log1p   WordCount.log1p  0.254319628               0 0.254319628
## NDSSName.my.fctr NDSSName.my.fctr  0.165445970               0 0.165445970
## WordCount.nexp     WordCount.nexp -0.053208396               0 0.053208396
## UniqueID                 UniqueID  0.011824920               1 0.011824920
## .rnorm                     .rnorm  0.008212201               0 0.008212201
```

```r
print(glb_feats_df <- orderBy(~-cor.y, 
          myfind_cor_features(feats_df=glb_feats_df, obs_df=glbObsTrn, rsp_var=glb_rsp_var,
                              nzv.freqCut=glb_nzv_freqCut, nzv.uniqueCut=glb_nzv_uniqueCut)))
```

```
## [1] "cor(WordCount.log1p, WordCount.root2)=0.8906"
## [1] "cor(Popular.fctr, WordCount.log1p)=0.2543"
## [1] "cor(Popular.fctr, WordCount.root2)=0.2921"
```

```
## Warning in myfind_cor_features(feats_df = glb_feats_df, obs_df =
## glbObsTrn, : Identified WordCount.log1p as highly correlated with
## WordCount.root2
```

```
##                                id        cor.y exclude.as.feat   cor.y.abs
## Popular                   Popular  1.000000000               1 1.000000000
## WordCount.root2   WordCount.root2  0.292120679               0 0.292120679
## WordCount               WordCount  0.257526549               1 0.257526549
## WordCount.log1p   WordCount.log1p  0.254319628               0 0.254319628
## NDSSName.my.fctr NDSSName.my.fctr  0.165445970               0 0.165445970
## UniqueID                 UniqueID  0.011824920               1 0.011824920
## .rnorm                     .rnorm  0.008212201               0 0.008212201
## WordCount.nexp     WordCount.nexp -0.053208396               0 0.053208396
##                       cor.high.X freqRatio percentUnique zeroVar   nzv
## Popular                     <NA>  4.976212    0.03061849   FALSE FALSE
## WordCount.root2             <NA>  2.315789   24.15799143   FALSE FALSE
## WordCount                   <NA>  2.315789   24.15799143   FALSE FALSE
## WordCount.log1p  WordCount.root2  2.315789   24.15799143   FALSE FALSE
## NDSSName.my.fctr            <NA>  1.348739    0.32149418   FALSE FALSE
## UniqueID                    <NA>  1.000000  100.00000000   FALSE FALSE
## .rnorm                      <NA>  1.000000  100.00000000   FALSE FALSE
## WordCount.nexp              <NA> 17.761364   11.32884262   FALSE FALSE
##                  is.cor.y.abs.low
## Popular                     FALSE
## WordCount.root2             FALSE
## WordCount                   FALSE
## WordCount.log1p             FALSE
## NDSSName.my.fctr            FALSE
## UniqueID                    FALSE
## .rnorm                      FALSE
## WordCount.nexp              FALSE
```

```r
plt_feats_df <- glb_feats_df
print(myplot_scatter(plt_feats_df, "percentUnique", "freqRatio", 
                     colorcol_name="nzv", jitter=TRUE) + 
          #geom_point(aes(shape=nzv)) +           
          geom_point() + 
          xlim(-5, 25) + 
          geom_hline(yintercept=glb_nzv_freqCut) +
          geom_vline(xintercept=glb_nzv_uniqueCut))
```

```
## Warning in myplot_scatter(plt_feats_df, "percentUnique", "freqRatio",
## colorcol_name = "nzv", : converting nzv to class:factor
```

```
## Warning: Removed 2 rows containing missing values (geom_point).
```

```
## Warning: Removed 2 rows containing missing values (geom_point).
```

```
## Warning: Removed 2 rows containing missing values (geom_point).
```

![](NYTBlogs3_base_files/figure-html/select.features-1.png) 

```r
print(subset(glb_feats_df, nzv))
```

```
##  [1] id               cor.y            exclude.as.feat  cor.y.abs       
##  [5] cor.high.X       freqRatio        percentUnique    zeroVar         
##  [9] nzv              is.cor.y.abs.low
## <0 rows> (or 0-length row.names)
```

```r
tmp_allobs_df <- 
    glbObsAll[, union(setdiff(names(glbObsAll), subset(glb_feats_df, nzv)$id),
                          glb_cluster_entropy_var)]
glbObsTrn <- subset(tmp_allobs_df, .src == "Train")
glbObsNew <- subset(tmp_allobs_df, .src == "Test")

glb_feats_df$interaction.feat <- NA
for (feat in names(glbFeatsInteractionOnly))
    glb_feats_df[glb_feats_df$id %in% feat, "interaction.feat"] <-
        glbFeatsInteractionOnly[[feat]]
        
#stop(here"); glb_to_sav(); glbObsAll <- sav_allobs_df
indep_vars <- subset(glb_feats_df, !nzv & (exclude.as.feat != 1))[, "id"]
numeric_indep_vars <- indep_vars[!grepl(".fctr", indep_vars, fixed=TRUE)]
glb_feats_df$shapiro.test.p.value <- NA
glb_feats_df[glb_feats_df$id %in% numeric_indep_vars, "shapiro.test.p.value"] <- 
    sapply(numeric_indep_vars, function(var) 
        shapiro.test(glbObsTrn[, var][1:min(5000, nrow(glbObsTrn))])$p.value)
not_nrml_feats_df <- glb_feats_df %>%
                        subset(!is.na(shapiro.test.p.value)) %>%
                    subset((shapiro.test.p.value < 0.05) || (id == ".rnorm")) %>%
                        arrange(shapiro.test.p.value)
row.names(not_nrml_feats_df) <- not_nrml_feats_df$id

#pltTrnObs <- glbObsTrn[, c("D.npnct05.log", ".rnorm")]
pltTrnObs <- glbObsTrn[, c(union(".rnorm",
                        not_nrml_feats_df$id[1:min(5, nrow(not_nrml_feats_df))]),
                            glb_cluster_entropy_var)]
gp <- myplot_violin(pltTrnObs, 
                    setdiff(names(pltTrnObs), glb_cluster_entropy_var), 
                    xcol_name = glb_cluster_entropy_var)
if ("variable" %in% names(pltTrnObs))
    gp <- gp + facet_wrap(~variable, scales = "free")
print(gp)
```

![](NYTBlogs3_base_files/figure-html/select.features-2.png) 

```r
#myplot_histogram(pltTrnObs, "D.npnct11.log", fill_col_name="sold", show_stats = TRUE)

myadjust_interaction_feats <- function(vars_vctr) {
    for (feat in subset(glb_feats_df, !is.na(interaction.feat))$id)
        if (feat %in% vars_vctr)
            vars_vctr <- union(setdiff(vars_vctr, feat), 
                paste0(glb_feats_df[glb_feats_df$id == feat, "interaction.feat"], ":",
                       feat))
    return(vars_vctr)
}

myrun_rfe <- function(obs_df, indep_vars, sizes = NULL) {
    rfe_obs_df <- myget_vectorized_obs_df(obs_df, glb_rsp_var, indep_vars)
    predictors_vctr <- setdiff(names(rfe_obs_df), glb_rsp_var)
    
    if (glb_is_regression)  rfeFuncs <- lmFuncs else {    
        rfeFuncs <- ldaFuncs
        
        # Delete non-variant columns
        predictors_unqLen <- sapply(predictors_vctr, function(col)
                                                length(unique(rfe_obs_df[, col])))
        predictors_vctr <- predictors_vctr[predictors_unqLen > 1]
        # Delete freqRatio >= 291
        #   plagiarized from caret:::nzv
        predictors_freqRatio <- 
            apply(rfe_obs_df[, predictors_vctr, FALSE], 2, function(data) {
                t <- table(data[!is.na(data)])
                if (length(t) <= 1) {
                    return(0)
                }
                w <- which.max(t)
                return(max(t, na.rm = TRUE)/max(t[-w], na.rm = TRUE))
            })
        predictors_vctr <- predictors_vctr[predictors_freqRatio < 172]
    }                        
    
    if (is.null(sizes))
        sizes <- tail(2 ^ (1:as.integer(log2(length(predictors_vctr)))), 5)
    
    rfe_control <- rfeControl(functions = rfeFuncs, method = "repeatedcv",
                             number = glb_rcv_n_folds, repeats = glb_rcv_n_repeats,
                              verbose = TRUE, returnResamp = "all",
        seeds = mygen_seeds(seeds_lst_len = 
                                (glb_rcv_n_folds * glb_rcv_n_repeats) + 1,
                            seeds_elmnt_lst_len = (length(sizes) + 1))
                            , allowParallel = FALSE
                            )
    set.seed(113)
    rfe_results <- rfe(rfe_obs_df[, predictors_vctr, FALSE], 
                       rfe_obs_df[, glb_rsp_var],
                       sizes = sizes, 
                       # metric = unlist(strsplit(glbMdlMetricsEval, "[.]"))[2],
#         maximize = ifelse(unlist(strsplit(glbMdlMetricsEval, "[.]"))[1] == "max",
#                                        TRUE, FALSE),
                       rfeControl = rfe_control)
    print(rfe_results)
    print(predictors(rfe_results))
    # print(plot(rfe_results, type=c("g", "o")))
    # print(plot(rfe_results))
    print(ggplot(rfe_results))

    return(rfe_results)
}

#stop(here"); glb_to_sav()
# shd .clusterid.fctr be excluded from this ? or include encoding of glb_category_var:.clusterid.fctr ?
indep_vars <- myadjust_interaction_feats(subset(glb_feats_df, 
                                                !nzv & (exclude.as.feat != 1))$id)
rfe_fit_results <- myrun_rfe(obs_df = glbObsFit, indep_vars = indep_vars, 
                             sizes = glbRFESizes[["RFE.X"]])
```

```
## +(rfe) fit Fold1.Rep1 size: 23 
## -(rfe) fit Fold1.Rep1 size: 23 
## +(rfe) imp Fold1.Rep1 
## -(rfe) imp Fold1.Rep1 
## +(rfe) fit Fold1.Rep1 size: 16 
## -(rfe) fit Fold1.Rep1 size: 16 
## +(rfe) fit Fold1.Rep1 size:  8 
## -(rfe) fit Fold1.Rep1 size:  8 
## +(rfe) fit Fold1.Rep1 size:  4 
## -(rfe) fit Fold1.Rep1 size:  4 
## +(rfe) fit Fold1.Rep1 size:  2 
## -(rfe) fit Fold1.Rep1 size:  2 
## +(rfe) fit Fold2.Rep1 size: 23 
## -(rfe) fit Fold2.Rep1 size: 23 
## +(rfe) imp Fold2.Rep1 
## -(rfe) imp Fold2.Rep1 
## +(rfe) fit Fold2.Rep1 size: 16 
## -(rfe) fit Fold2.Rep1 size: 16 
## +(rfe) fit Fold2.Rep1 size:  8 
## -(rfe) fit Fold2.Rep1 size:  8 
## +(rfe) fit Fold2.Rep1 size:  4 
## -(rfe) fit Fold2.Rep1 size:  4 
## +(rfe) fit Fold2.Rep1 size:  2 
## -(rfe) fit Fold2.Rep1 size:  2 
## +(rfe) fit Fold3.Rep1 size: 23 
## -(rfe) fit Fold3.Rep1 size: 23 
## +(rfe) imp Fold3.Rep1 
## -(rfe) imp Fold3.Rep1 
## +(rfe) fit Fold3.Rep1 size: 16 
## -(rfe) fit Fold3.Rep1 size: 16 
## +(rfe) fit Fold3.Rep1 size:  8 
## -(rfe) fit Fold3.Rep1 size:  8 
## +(rfe) fit Fold3.Rep1 size:  4 
## -(rfe) fit Fold3.Rep1 size:  4 
## +(rfe) fit Fold3.Rep1 size:  2 
## -(rfe) fit Fold3.Rep1 size:  2 
## +(rfe) fit Fold1.Rep2 size: 23 
## -(rfe) fit Fold1.Rep2 size: 23 
## +(rfe) imp Fold1.Rep2 
## -(rfe) imp Fold1.Rep2 
## +(rfe) fit Fold1.Rep2 size: 16 
## -(rfe) fit Fold1.Rep2 size: 16 
## +(rfe) fit Fold1.Rep2 size:  8 
## -(rfe) fit Fold1.Rep2 size:  8 
## +(rfe) fit Fold1.Rep2 size:  4 
## -(rfe) fit Fold1.Rep2 size:  4 
## +(rfe) fit Fold1.Rep2 size:  2 
## -(rfe) fit Fold1.Rep2 size:  2 
## +(rfe) fit Fold2.Rep2 size: 23 
## -(rfe) fit Fold2.Rep2 size: 23 
## +(rfe) imp Fold2.Rep2 
## -(rfe) imp Fold2.Rep2 
## +(rfe) fit Fold2.Rep2 size: 16 
## -(rfe) fit Fold2.Rep2 size: 16 
## +(rfe) fit Fold2.Rep2 size:  8 
## -(rfe) fit Fold2.Rep2 size:  8 
## +(rfe) fit Fold2.Rep2 size:  4 
## -(rfe) fit Fold2.Rep2 size:  4 
## +(rfe) fit Fold2.Rep2 size:  2 
## -(rfe) fit Fold2.Rep2 size:  2 
## +(rfe) fit Fold3.Rep2 size: 23 
## -(rfe) fit Fold3.Rep2 size: 23 
## +(rfe) imp Fold3.Rep2 
## -(rfe) imp Fold3.Rep2 
## +(rfe) fit Fold3.Rep2 size: 16 
## -(rfe) fit Fold3.Rep2 size: 16 
## +(rfe) fit Fold3.Rep2 size:  8 
## -(rfe) fit Fold3.Rep2 size:  8 
## +(rfe) fit Fold3.Rep2 size:  4 
## -(rfe) fit Fold3.Rep2 size:  4 
## +(rfe) fit Fold3.Rep2 size:  2 
## -(rfe) fit Fold3.Rep2 size:  2 
## +(rfe) fit Fold1.Rep3 size: 23 
## -(rfe) fit Fold1.Rep3 size: 23 
## +(rfe) imp Fold1.Rep3 
## -(rfe) imp Fold1.Rep3 
## +(rfe) fit Fold1.Rep3 size: 16 
## -(rfe) fit Fold1.Rep3 size: 16 
## +(rfe) fit Fold1.Rep3 size:  8 
## -(rfe) fit Fold1.Rep3 size:  8 
## +(rfe) fit Fold1.Rep3 size:  4 
## -(rfe) fit Fold1.Rep3 size:  4 
## +(rfe) fit Fold1.Rep3 size:  2 
## -(rfe) fit Fold1.Rep3 size:  2 
## +(rfe) fit Fold2.Rep3 size: 23 
## -(rfe) fit Fold2.Rep3 size: 23 
## +(rfe) imp Fold2.Rep3 
## -(rfe) imp Fold2.Rep3 
## +(rfe) fit Fold2.Rep3 size: 16 
## -(rfe) fit Fold2.Rep3 size: 16 
## +(rfe) fit Fold2.Rep3 size:  8 
## -(rfe) fit Fold2.Rep3 size:  8 
## +(rfe) fit Fold2.Rep3 size:  4 
## -(rfe) fit Fold2.Rep3 size:  4 
## +(rfe) fit Fold2.Rep3 size:  2 
## -(rfe) fit Fold2.Rep3 size:  2 
## +(rfe) fit Fold3.Rep3 size: 23 
## -(rfe) fit Fold3.Rep3 size: 23 
## +(rfe) imp Fold3.Rep3 
## -(rfe) imp Fold3.Rep3 
## +(rfe) fit Fold3.Rep3 size: 16 
## -(rfe) fit Fold3.Rep3 size: 16 
## +(rfe) fit Fold3.Rep3 size:  8 
## -(rfe) fit Fold3.Rep3 size:  8 
## +(rfe) fit Fold3.Rep3 size:  4 
## -(rfe) fit Fold3.Rep3 size:  4 
## +(rfe) fit Fold3.Rep3 size:  2 
## -(rfe) fit Fold3.Rep3 size:  2 
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
##          8   0.9304 0.75904   0.004609 0.01699         
##         16   0.9301 0.75823   0.004601 0.01698         
##         23   0.9325 0.76819   0.004975 0.01772        *
## 
## The top 5 variables (out of 23):
##    WordCount.log1p, WordCount.root2, WordCount.nexp, NDSSName.my.fctrOpEd#Opinion#, NDSSName.my.fctrScience#Health#
## 
##  [1] "WordCount.log1p"                                   
##  [2] "WordCount.root2"                                   
##  [3] "WordCount.nexp"                                    
##  [4] "NDSSName.my.fctrOpEd#Opinion#"                     
##  [5] "NDSSName.my.fctrScience#Health#"                   
##  [6] "NDSSName.my.fctrBusiness#Crosswords/Games#"        
##  [7] "NDSSName.my.fctrStyles#U.S.#"                      
##  [8] ".rnorm"                                            
##  [9] "NDSSName.my.fctrmyOther"                           
## [10] "NDSSName.my.fctr#Opinion#RoomForDebate"            
## [11] "NDSSName.my.fctrBusiness#Technology#"              
## [12] "NDSSName.my.fctrMetro#N.Y./Region#"                
## [13] "NDSSName.my.fctrTravel#Travel#"                    
## [14] "NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness"
## [15] "NDSSName.my.fctr#Multimedia#"                      
## [16] "NDSSName.my.fctrStyles##Fashion"                   
## [17] "NDSSName.my.fctrForeign#World#"                    
## [18] "NDSSName.my.fctrForeign#World#AsiaPacific"         
## [19] "NDSSName.my.fctr#U.S.#Education"                   
## [20] "NDSSName.my.fctrCulture#Arts#"                     
## [21] "NDSSName.my.fctrBusiness#BusinessDay#Dealbook"     
## [22] "NDSSName.my.fctr##"                                
## [23] "NDSSName.my.fctrTStyle##"
```

![](NYTBlogs3_base_files/figure-html/select.features-3.png) 

```r
# print(all.equal(rfe_results[-which(names(rfe_results) == "times")], 
#                 sav_rfe_results[-which(names(sav_rfe_results) == "times")]))

# require(mRMRe)
# indep_vars_vctr <- subset(glb_feats_df, !nzv &
#                                         (exclude.as.feat != 1))[, "id"]
# indep_vars_vctr <- setdiff(indep_vars_vctr, 
#                     myfind_fctr_cols_df(glbObsTrn[, c(glb_rsp_var, indep_vars_vctr)]))
# tmp_trnobs_df <- glbObsTrn[, c(glb_rsp_var, indep_vars_vctr)]
# tmp_trnobs_df$biddable <- as.numeric(tmp_trnobs_df$biddable)
# dd <- mRMR.data(data = tmp_trnobs_df)
# mRMRe.fltr <- mRMR.classic(data = dd, target_indices = c(1), feature_count = 10)
# print(solutions(mRMRe.fltr)[[1]])
# print(apply(solutions(mRMRe.fltr)[[1]], 2, function(x, y) { return(y[x]) },
#             y=featureNames(dd)))
# print(featureNames(dd)[solutions(mRMRe.fltr)[[1]]])
# print(mRMRe.fltr@filters); print(mRMRe.fltr@scores)

mycheck_problem_data(glbObsAll, featsExclude = glbFeatsExclude, 
                     fctrMaxUniqVals = glbFctrMaxUniqVals, terminate = TRUE)
```

```
## [1] "numeric data missing in glbObsAll: "
##      Popular Popular.fctr 
##         1870         1870 
## [1] "numeric data w/ 0s in glbObsAll: "
##       WordCount         Popular WordCount.log1p WordCount.root2 
##             109            5439             109             109 
##  WordCount.nexp 
##            2044 
## [1] "numeric data w/ Infs in glbObsAll: "
## named integer(0)
## [1] "numeric data w/ NaNs in glbObsAll: "
## named integer(0)
## [1] "string data missing in glbObsAll: "
##       NewsDesk    SectionName SubsectionName       Headline        Snippet 
##           2408           2899           6176              0             13 
##       Abstract        PubDate    NDSSName.my           .lcn 
##             17              0              0           1870
```

```r
# glbObsAll %>% filter(is.na(Married.fctr)) %>% tbl_df()
# glbObsAll %>% count(Married.fctr)
# levels(glbObsAll$Married.fctr)

print("glb_feats_df:");   print(dim(glb_feats_df))
```

```
## [1] "glb_feats_df:"
```

```
## [1]  8 12
```

```r
sav_feats_df <- glb_feats_df
glb_feats_df <- sav_feats_df

glb_feats_df[, "rsp_var_raw"] <- FALSE
glb_feats_df[glb_feats_df$id == glb_rsp_var_raw, "rsp_var_raw"] <- TRUE 
glb_feats_df$exclude.as.feat <- (glb_feats_df$exclude.as.feat == 1)
if (!is.null(glb_id_var) && glb_id_var != ".rownames")
    glb_feats_df[glb_feats_df$id %in% glb_id_var, "id_var"] <- TRUE 
add_feats_df <- data.frame(id=glb_rsp_var, exclude.as.feat=TRUE, rsp_var=TRUE)
row.names(add_feats_df) <- add_feats_df$id; print(add_feats_df)
```

```
##                        id exclude.as.feat rsp_var
## Popular.fctr Popular.fctr            TRUE    TRUE
```

```r
glb_feats_df <- myrbind_df(glb_feats_df, add_feats_df)
if (glb_id_var != ".rownames")
    print(subset(glb_feats_df, rsp_var_raw | rsp_var | id_var)) else
    print(subset(glb_feats_df, rsp_var_raw | rsp_var))    
```

```
##                        id      cor.y exclude.as.feat  cor.y.abs cor.high.X
## Popular           Popular 1.00000000            TRUE 1.00000000       <NA>
## UniqueID         UniqueID 0.01182492            TRUE 0.01182492       <NA>
## Popular.fctr Popular.fctr         NA            TRUE         NA       <NA>
##              freqRatio percentUnique zeroVar   nzv is.cor.y.abs.low
## Popular       4.976212    0.03061849   FALSE FALSE            FALSE
## UniqueID      1.000000  100.00000000   FALSE FALSE            FALSE
## Popular.fctr        NA            NA      NA    NA               NA
##              interaction.feat shapiro.test.p.value rsp_var_raw id_var
## Popular                    NA                   NA        TRUE     NA
## UniqueID                   NA                   NA       FALSE   TRUE
## Popular.fctr               NA                   NA          NA     NA
##              rsp_var
## Popular           NA
## UniqueID          NA
## Popular.fctr    TRUE
```

```r
print("glb_feats_df vs. glbObsAll: "); 
```

```
## [1] "glb_feats_df vs. glbObsAll: "
```

```r
print(setdiff(glb_feats_df$id, names(glbObsAll)))
```

```
## character(0)
```

```r
print("glbObsAll vs. glb_feats_df: "); 
```

```
## [1] "glbObsAll vs. glb_feats_df: "
```

```r
# Ensure these are only chr vars
print(setdiff(setdiff(names(glbObsAll), glb_feats_df$id), 
                myfind_chr_cols_df(glbObsAll)))
```

```
## character(0)
```

```r
if (glb_save_envir)
    save(glb_feats_df, 
         glbObsAll, #glbObsTrn, glbObsFit, glbObsOOB, glbObsNew,
         file=paste0(glb_out_pfx, "selfts_dsk.RData"))
# load(paste0(glb_out_pfx, "blddfs_dsk.RData"))

# if (!all.equal(tmp_feats_df, glb_feats_df))
#     stop("glb_feats_df r/w not working")
# if (!all.equal(tmp_entity_df, glbObsAll))
#     stop("glbObsAll r/w not working")

glb_chunks_df <- myadd_chunk(glb_chunks_df, "fit.models", major.inc=TRUE)
```

```
##              label step_major step_minor label_minor    bgn    end elapsed
## 9  select.features          5          0           0 25.428 35.237   9.809
## 10      fit.models          6          0           0 35.238     NA      NA
```

## Step `6.0: fit models`

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
if (!is.null(glb_Baseline_mdl_var)) 
    ret_lst <- myfit_mdl(mdl_id="Baseline", 
                         model_method="mybaseln_classfr",
                        indep_vars_vctr=glb_Baseline_mdl_var,
                        rsp_var=glb_rsp_var, rsp_var_out=glb_rsp_var_out,
                        fit_df=glbObsFit, OOB_df=glbObsOOB)

# Most Frequent Outcome "MFO" model: mean(y) for regression
#   Not using caret's nullModel since model stats not avl
#   Cannot use rpart for multinomial classification since it predicts non-MFO
ret_lst <- myfit_mdl(mdl_specs_lst = myinit_mdl_specs_lst(mdl_specs_lst = list(
    id.prefix = "MFO", type = glb_model_type, trainControl.method = "none",
    train.method = ifelse(glb_is_regression, "lm", "myMFO_classfr"))),
                        indep_vars = ".rnorm", rsp_var = glb_rsp_var,
                        fit_df = glbObsFit, OOB_df = glbObsOOB)
```

```
## [1] "fitting model: MFO.myMFO_classfr"
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

![](NYTBlogs3_base_files/figure-html/fit.models_0-1.png) 

```
##   Popular.fctr Popular.fctr.predict.MFO.myMFO_classfr.Y
## 1            N                                     3941
## 2            Y                                      863
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

![](NYTBlogs3_base_files/figure-html/fit.models_0-2.png) 

```
##   Popular.fctr Popular.fctr.predict.MFO.myMFO_classfr.Y
## 1            N                                     1498
## 2            Y                                      230
##          Prediction
## Reference    N    Y
##         N    0 1498
##         Y    0  230
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##      0.1331019      0.0000000      0.1174298      0.1500310      0.8668981 
## AccuracyPValue  McnemarPValue 
##      1.0000000      0.0000000 
##                  id  feats max.nTuningRuns min.elapsedtime.everything
## 1 MFO.myMFO_classfr .rnorm               0                      0.269
##   min.elapsedtime.final max.AUCpROC.fit max.Sens.fit max.Spec.fit
## 1                 0.003             0.5            1            0
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
if (glb_is_classification)
    # "random" model - only for classification; 
    #   none needed for regression since it is same as MFO
    ret_lst <- myfit_mdl(mdl_specs_lst = myinit_mdl_specs_lst(mdl_specs_lst = list(
        id.prefix = "Random", type = glb_model_type, trainControl.method = "none",
        train.method = "myrandom_classfr")),
                        indep_vars = ".rnorm", rsp_var = glb_rsp_var,
                        fit_df = glbObsFit, OOB_df = glbObsOOB)
```

```
## [1] "fitting model: Random.myrandom_classfr"
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

![](NYTBlogs3_base_files/figure-html/fit.models_0-3.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-4.png) 

```
##   Popular.fctr Popular.fctr.predict.Random.myrandom_classfr.Y
## 1            N                                           3941
## 2            Y                                            863
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

![](NYTBlogs3_base_files/figure-html/fit.models_0-5.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-6.png) 

```
##   Popular.fctr Popular.fctr.predict.Random.myrandom_classfr.Y
## 1            N                                           1498
## 2            Y                                            230
##          Prediction
## Reference    N    Y
##         N    0 1498
##         Y    0  230
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##      0.1331019      0.0000000      0.1174298      0.1500310      0.8668981 
## AccuracyPValue  McnemarPValue 
##      1.0000000      0.0000000 
##                        id  feats max.nTuningRuns
## 1 Random.myrandom_classfr .rnorm               0
##   min.elapsedtime.everything min.elapsedtime.final max.AUCpROC.fit
## 1                      0.269                 0.001       0.4990604
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
#     ret_lst <- myfit_mdl(mdl_id = "Random", model_method = "myrandom_classfr",
#                             model_type = glb_model_type,                         
#                             indep_vars_vctr = ".rnorm",
#                             rsp_var = glb_rsp_var, rsp_var_out = glb_rsp_var_out,
#                             fit_df = glbObsFit, OOB_df = glbObsOOB)

# Max.cor.Y
#   Check impact of cv
#       rpart is not a good candidate since caret does not optimize cp (only tuning parameter of rpart) well
ret_lst <- myfit_mdl(mdl_specs_lst=myinit_mdl_specs_lst(mdl_specs_lst=list(
    id.prefix="Max.cor.Y.rcv.1X1", type=glb_model_type, trainControl.method="none",
    train.method="glmnet")),
                    indep_vars=max_cor_y_x_vars, rsp_var=glb_rsp_var, 
                    fit_df=glbObsFit, OOB_df=glbObsOOB)
```

```
## [1] "fitting model: Max.cor.Y.rcv.1X1.glmnet"
## [1] "    indep_vars: WordCount.root2,NDSSName.my.fctr"
```

```
## Loading required package: glmnet
## Loading required package: Matrix
## Loaded glmnet 2.0-2
```

```
## Fitting alpha = 0.1, lambda = 0.00434 on full training set
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-7.png) 

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
##                                        (Intercept) 
##                                        -4.57159198 
##                       NDSSName.my.fctr#Multimedia# 
##                                        -1.22219085 
##             NDSSName.my.fctr#Opinion#RoomForDebate 
##                                        -3.46072453 
##           NDSSName.my.fctr#Opinion#ThePublicEditor 
##                                         4.06871185 
##                    NDSSName.my.fctr#U.S.#Education 
##                                        -1.89443632 
##      NDSSName.my.fctrBusiness#BusinessDay#Dealbook 
##                                        -0.22472818 
## NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness 
##                                        -0.95537118 
##         NDSSName.my.fctrBusiness#Crosswords/Games# 
##                                         4.55408513 
##               NDSSName.my.fctrBusiness#Technology# 
##                                         0.77368538 
##                      NDSSName.my.fctrCulture#Arts# 
##                                        -0.09465691 
##                     NDSSName.my.fctrForeign#World# 
##                                        -1.45528874 
##          NDSSName.my.fctrForeign#World#AsiaPacific 
##                                        -1.60117505 
##                 NDSSName.my.fctrMetro#N.Y./Region# 
##                                         0.01563989 
##                      NDSSName.my.fctrOpEd#Opinion# 
##                                         4.51696382 
##                    NDSSName.my.fctrScience#Health# 
##                                         3.51595317 
##                    NDSSName.my.fctrStyles##Fashion 
##                                        -1.85948925 
##                       NDSSName.my.fctrStyles#U.S.# 
##                                         3.27995325 
##                           NDSSName.my.fctrTStyle## 
##                                        -1.54110404 
##                     NDSSName.my.fctrTravel#Travel# 
##                                        -1.41940605 
##                            NDSSName.my.fctrmyOther 
##                                        -1.90156922 
##                                    WordCount.root2 
##                                         0.08434378 
## [1] "max lambda < lambdaOpt:"
##                                        (Intercept) 
##                                        -4.60394059 
##                       NDSSName.my.fctr#Multimedia# 
##                                        -1.25163328 
##             NDSSName.my.fctr#Opinion#RoomForDebate 
##                                        -3.55521332 
##           NDSSName.my.fctr#Opinion#ThePublicEditor 
##                                         4.09217313 
##                    NDSSName.my.fctr#U.S.#Education 
##                                        -1.96172971 
##      NDSSName.my.fctrBusiness#BusinessDay#Dealbook 
##                                        -0.22495986 
## NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness 
##                                        -0.96836050 
##         NDSSName.my.fctrBusiness#Crosswords/Games# 
##                                         4.58120497 
##               NDSSName.my.fctrBusiness#Technology# 
##                                         0.78504703 
##                      NDSSName.my.fctrCulture#Arts# 
##                                        -0.09069661 
##                     NDSSName.my.fctrForeign#World# 
##                                        -1.51061232 
##          NDSSName.my.fctrForeign#World#AsiaPacific 
##                                        -1.63313235 
##                 NDSSName.my.fctrMetro#N.Y./Region# 
##                                         0.02466697 
##                      NDSSName.my.fctrOpEd#Opinion# 
##                                         4.54361134 
##                    NDSSName.my.fctrScience#Health# 
##                                         3.53210055 
##                    NDSSName.my.fctrStyles##Fashion 
##                                        -1.92188290 
##                       NDSSName.my.fctrStyles#U.S.# 
##                                         3.29488750 
##                           NDSSName.my.fctrTStyle## 
##                                        -1.57788931 
##                     NDSSName.my.fctrTravel#Travel# 
##                                        -1.47368131 
##                            NDSSName.my.fctrmyOther 
##                                        -1.97357582 
##                                    WordCount.root2 
##                                         0.08537319
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-8.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-9.png) 

```
##   Popular.fctr Popular.fctr.predict.Max.cor.Y.rcv.1X1.glmnet.N
## 1            N                                            3796
## 2            Y                                             177
##   Popular.fctr.predict.Max.cor.Y.rcv.1X1.glmnet.Y
## 1                                             145
## 2                                             686
##          Prediction
## Reference    N    Y
##         N 3796  145
##         Y  177  686
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   9.329725e-01   7.692476e-01   9.255302e-01   9.398832e-01   8.203580e-01 
## AccuracyPValue  McnemarPValue 
##  1.026390e-114   8.406670e-02
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-10.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-11.png) 

```
##   Popular.fctr Popular.fctr.predict.Max.cor.Y.rcv.1X1.glmnet.N
## 1            N                                            1151
## 2            Y                                              67
##   Popular.fctr.predict.Max.cor.Y.rcv.1X1.glmnet.Y
## 1                                             347
## 2                                             163
##          Prediction
## Reference    N    Y
##         N 1151  347
##         Y   67  163
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   7.604167e-01   3.148374e-01   7.395703e-01   7.803749e-01   8.668981e-01 
## AccuracyPValue  McnemarPValue 
##   1.000000e+00   8.593187e-43 
##                         id                            feats
## 1 Max.cor.Y.rcv.1X1.glmnet WordCount.root2,NDSSName.my.fctr
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1               0                      0.972                 0.273
##   max.AUCpROC.fit max.Sens.fit max.Spec.fit max.AUCROCR.fit
## 1       0.8790544    0.9632073    0.7949015       0.9608594
##   opt.prob.threshold.fit max.f.score.fit max.Accuracy.fit
## 1                    0.5       0.8099174        0.9329725
##   max.AccuracyLower.fit max.AccuracyUpper.fit max.Kappa.fit
## 1             0.9255302             0.9398832     0.7692476
##   max.AUCpROC.OOB max.Sens.OOB max.Spec.OOB max.AUCROCR.OOB
## 1       0.5962443    0.9098798    0.2826087       0.8116126
##   opt.prob.threshold.OOB max.f.score.OOB max.Accuracy.OOB
## 1                    0.1       0.4405405        0.7604167
##   max.AccuracyLower.OOB max.AccuracyUpper.OOB max.Kappa.OOB
## 1             0.7395703             0.7803749     0.3148374
```

```r
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
```

```
## [1] "fitting model: Max.cor.Y.rcv.3X1.glmnet"
## [1] "    indep_vars: WordCount.root2,NDSSName.my.fctr"
## Aggregating results
## Selecting tuning parameters
## Fitting alpha = 0.325, lambda = 0.0201 on full training set
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-12.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-13.png) 

```
##             Length Class      Mode     
## a0            99   -none-     numeric  
## beta        2079   dgCMatrix  S4       
## df            99   -none-     numeric  
## dim            2   -none-     numeric  
## lambda        99   -none-     numeric  
## dev.ratio     99   -none-     numeric  
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
##                                        (Intercept) 
##                                        -3.89350373 
##                       NDSSName.my.fctr#Multimedia# 
##                                        -0.01916344 
##             NDSSName.my.fctr#Opinion#RoomForDebate 
##                                        -1.18453357 
##           NDSSName.my.fctr#Opinion#ThePublicEditor 
##                                         3.21701058 
##                    NDSSName.my.fctr#U.S.#Education 
##                                        -0.47679040 
## NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness 
##                                        -0.09891374 
##         NDSSName.my.fctrBusiness#Crosswords/Games# 
##                                         3.87281404 
##               NDSSName.my.fctrBusiness#Technology# 
##                                         0.40965256 
##                     NDSSName.my.fctrForeign#World# 
##                                        -0.05114617 
##          NDSSName.my.fctrForeign#World#AsiaPacific 
##                                        -0.47464340 
##                      NDSSName.my.fctrOpEd#Opinion# 
##                                         3.95214357 
##                    NDSSName.my.fctrScience#Health# 
##                                         3.14232408 
##                    NDSSName.my.fctrStyles##Fashion 
##                                        -0.31867093 
##                       NDSSName.my.fctrStyles#U.S.# 
##                                         2.92610567 
##                           NDSSName.my.fctrTStyle## 
##                                        -0.60538025 
##                                    WordCount.root2 
##                                         0.05783392 
## [1] "max lambda < lambdaOpt:"
##                                        (Intercept) 
##                                        -3.95644632 
##                       NDSSName.my.fctr#Multimedia# 
##                                        -0.07182859 
##             NDSSName.my.fctr#Opinion#RoomForDebate 
##                                        -1.30034382 
##           NDSSName.my.fctr#Opinion#ThePublicEditor 
##                                         3.29694236 
##                    NDSSName.my.fctr#U.S.#Education 
##                                        -0.53415905 
## NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness 
##                                        -0.14259759 
##         NDSSName.my.fctrBusiness#Crosswords/Games# 
##                                         3.94231812 
##               NDSSName.my.fctrBusiness#Technology# 
##                                         0.45657914 
##                     NDSSName.my.fctrForeign#World# 
##                                        -0.10021084 
##          NDSSName.my.fctrForeign#World#AsiaPacific 
##                                        -0.53077048 
##                      NDSSName.my.fctrOpEd#Opinion# 
##                                         4.01324666 
##                    NDSSName.my.fctrScience#Health# 
##                                         3.18936803 
##                    NDSSName.my.fctrStyles##Fashion 
##                                        -0.38069674 
##                       NDSSName.my.fctrStyles#U.S.# 
##                                         2.97176051 
##                           NDSSName.my.fctrTStyle## 
##                                        -0.64837249 
##                                    WordCount.root2 
##                                         0.05978318
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-14.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-15.png) 

```
##   Popular.fctr Popular.fctr.predict.Max.cor.Y.rcv.3X1.glmnet.N
## 1            N                                            3796
## 2            Y                                             177
##   Popular.fctr.predict.Max.cor.Y.rcv.3X1.glmnet.Y
## 1                                             145
## 2                                             686
##          Prediction
## Reference    N    Y
##         N 3796  145
##         Y  177  686
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   9.329725e-01   7.692476e-01   9.255302e-01   9.398832e-01   8.203580e-01 
## AccuracyPValue  McnemarPValue 
##  1.026390e-114   8.406670e-02
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-16.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-17.png) 

```
##   Popular.fctr Popular.fctr.predict.Max.cor.Y.rcv.3X1.glmnet.N
## 1            N                                            1146
## 2            Y                                              67
##   Popular.fctr.predict.Max.cor.Y.rcv.3X1.glmnet.Y
## 1                                             352
## 2                                             163
##          Prediction
## Reference    N    Y
##         N 1146  352
##         Y   67  163
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   7.575231e-01   3.107477e-01   7.365992e-01   7.775689e-01   8.668981e-01 
## AccuracyPValue  McnemarPValue 
##   1.000000e+00   9.066396e-44 
##                         id                            feats
## 1 Max.cor.Y.rcv.3X1.glmnet WordCount.root2,NDSSName.my.fctr
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1              25                      2.826                  0.27
##   max.AUCpROC.fit max.Sens.fit max.Spec.fit max.AUCROCR.fit
## 1       0.8767919     0.964476    0.7891078       0.9582555
##   opt.prob.threshold.fit max.f.score.fit max.Accuracy.fit
## 1                    0.4       0.8099174        0.9335973
##   max.AccuracyLower.fit max.AccuracyUpper.fit max.Kappa.fit
## 1             0.9255302             0.9398832     0.7691678
##   max.AUCpROC.OOB max.Sens.OOB max.Spec.OOB max.AUCROCR.OOB
## 1       0.5962443    0.9098798    0.2826087       0.8067975
##   opt.prob.threshold.OOB max.f.score.OOB max.Accuracy.OOB
## 1                    0.1       0.4375839        0.7575231
##   max.AccuracyLower.OOB max.AccuracyUpper.OOB max.Kappa.OOB
## 1             0.7365992             0.7775689     0.3107477
##   max.AccuracySD.fit max.KappaSD.fit
## 1        0.007015493      0.02403706
## [1] "fitting model: Max.cor.Y.rcv.3X3.glmnet"
## [1] "    indep_vars: WordCount.root2,NDSSName.my.fctr"
## Aggregating results
## Selecting tuning parameters
## Fitting alpha = 0.325, lambda = 0.0201 on full training set
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-18.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-19.png) 

```
##             Length Class      Mode     
## a0            99   -none-     numeric  
## beta        2079   dgCMatrix  S4       
## df            99   -none-     numeric  
## dim            2   -none-     numeric  
## lambda        99   -none-     numeric  
## dev.ratio     99   -none-     numeric  
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
##                                        (Intercept) 
##                                        -3.89350373 
##                       NDSSName.my.fctr#Multimedia# 
##                                        -0.01916344 
##             NDSSName.my.fctr#Opinion#RoomForDebate 
##                                        -1.18453357 
##           NDSSName.my.fctr#Opinion#ThePublicEditor 
##                                         3.21701058 
##                    NDSSName.my.fctr#U.S.#Education 
##                                        -0.47679040 
## NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness 
##                                        -0.09891374 
##         NDSSName.my.fctrBusiness#Crosswords/Games# 
##                                         3.87281404 
##               NDSSName.my.fctrBusiness#Technology# 
##                                         0.40965256 
##                     NDSSName.my.fctrForeign#World# 
##                                        -0.05114617 
##          NDSSName.my.fctrForeign#World#AsiaPacific 
##                                        -0.47464340 
##                      NDSSName.my.fctrOpEd#Opinion# 
##                                         3.95214357 
##                    NDSSName.my.fctrScience#Health# 
##                                         3.14232408 
##                    NDSSName.my.fctrStyles##Fashion 
##                                        -0.31867093 
##                       NDSSName.my.fctrStyles#U.S.# 
##                                         2.92610567 
##                           NDSSName.my.fctrTStyle## 
##                                        -0.60538025 
##                                    WordCount.root2 
##                                         0.05783392 
## [1] "max lambda < lambdaOpt:"
##                                        (Intercept) 
##                                        -3.95644632 
##                       NDSSName.my.fctr#Multimedia# 
##                                        -0.07182859 
##             NDSSName.my.fctr#Opinion#RoomForDebate 
##                                        -1.30034382 
##           NDSSName.my.fctr#Opinion#ThePublicEditor 
##                                         3.29694236 
##                    NDSSName.my.fctr#U.S.#Education 
##                                        -0.53415905 
## NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness 
##                                        -0.14259759 
##         NDSSName.my.fctrBusiness#Crosswords/Games# 
##                                         3.94231812 
##               NDSSName.my.fctrBusiness#Technology# 
##                                         0.45657914 
##                     NDSSName.my.fctrForeign#World# 
##                                        -0.10021084 
##          NDSSName.my.fctrForeign#World#AsiaPacific 
##                                        -0.53077048 
##                      NDSSName.my.fctrOpEd#Opinion# 
##                                         4.01324666 
##                    NDSSName.my.fctrScience#Health# 
##                                         3.18936803 
##                    NDSSName.my.fctrStyles##Fashion 
##                                        -0.38069674 
##                       NDSSName.my.fctrStyles#U.S.# 
##                                         2.97176051 
##                           NDSSName.my.fctrTStyle## 
##                                        -0.64837249 
##                                    WordCount.root2 
##                                         0.05978318
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-20.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-21.png) 

```
##   Popular.fctr Popular.fctr.predict.Max.cor.Y.rcv.3X3.glmnet.N
## 1            N                                            3796
## 2            Y                                             177
##   Popular.fctr.predict.Max.cor.Y.rcv.3X3.glmnet.Y
## 1                                             145
## 2                                             686
##          Prediction
## Reference    N    Y
##         N 3796  145
##         Y  177  686
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   9.329725e-01   7.692476e-01   9.255302e-01   9.398832e-01   8.203580e-01 
## AccuracyPValue  McnemarPValue 
##  1.026390e-114   8.406670e-02
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-22.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-23.png) 

```
##   Popular.fctr Popular.fctr.predict.Max.cor.Y.rcv.3X3.glmnet.N
## 1            N                                            1146
## 2            Y                                              67
##   Popular.fctr.predict.Max.cor.Y.rcv.3X3.glmnet.Y
## 1                                             352
## 2                                             163
##          Prediction
## Reference    N    Y
##         N 1146  352
##         Y   67  163
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   7.575231e-01   3.107477e-01   7.365992e-01   7.775689e-01   8.668981e-01 
## AccuracyPValue  McnemarPValue 
##   1.000000e+00   9.066396e-44 
##                         id                            feats
## 1 Max.cor.Y.rcv.3X3.glmnet WordCount.root2,NDSSName.my.fctr
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1              25                      4.817                  0.27
##   max.AUCpROC.fit max.Sens.fit max.Spec.fit max.AUCROCR.fit
## 1       0.8767919     0.964476    0.7891078       0.9582555
##   opt.prob.threshold.fit max.f.score.fit max.Accuracy.fit
## 1                    0.4       0.8099174        0.9333193
##   max.AccuracyLower.fit max.AccuracyUpper.fit max.Kappa.fit
## 1             0.9255302             0.9398832     0.7690803
##   max.AUCpROC.OOB max.Sens.OOB max.Spec.OOB max.AUCROCR.OOB
## 1       0.5962443    0.9098798    0.2826087       0.8067975
##   opt.prob.threshold.OOB max.f.score.OOB max.Accuracy.OOB
## 1                    0.1       0.4375839        0.7575231
##   max.AccuracyLower.OOB max.AccuracyUpper.OOB max.Kappa.OOB
## 1             0.7365992             0.7775689     0.3107477
##   max.AccuracySD.fit max.KappaSD.fit
## 1        0.005178375      0.01754365
## [1] "fitting model: Max.cor.Y.rcv.3X5.glmnet"
## [1] "    indep_vars: WordCount.root2,NDSSName.my.fctr"
## Aggregating results
## Selecting tuning parameters
## Fitting alpha = 0.325, lambda = 0.0201 on full training set
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-24.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-25.png) 

```
##             Length Class      Mode     
## a0            99   -none-     numeric  
## beta        2079   dgCMatrix  S4       
## df            99   -none-     numeric  
## dim            2   -none-     numeric  
## lambda        99   -none-     numeric  
## dev.ratio     99   -none-     numeric  
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
##                                        (Intercept) 
##                                        -3.89350373 
##                       NDSSName.my.fctr#Multimedia# 
##                                        -0.01916344 
##             NDSSName.my.fctr#Opinion#RoomForDebate 
##                                        -1.18453357 
##           NDSSName.my.fctr#Opinion#ThePublicEditor 
##                                         3.21701058 
##                    NDSSName.my.fctr#U.S.#Education 
##                                        -0.47679040 
## NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness 
##                                        -0.09891374 
##         NDSSName.my.fctrBusiness#Crosswords/Games# 
##                                         3.87281404 
##               NDSSName.my.fctrBusiness#Technology# 
##                                         0.40965256 
##                     NDSSName.my.fctrForeign#World# 
##                                        -0.05114617 
##          NDSSName.my.fctrForeign#World#AsiaPacific 
##                                        -0.47464340 
##                      NDSSName.my.fctrOpEd#Opinion# 
##                                         3.95214357 
##                    NDSSName.my.fctrScience#Health# 
##                                         3.14232408 
##                    NDSSName.my.fctrStyles##Fashion 
##                                        -0.31867093 
##                       NDSSName.my.fctrStyles#U.S.# 
##                                         2.92610567 
##                           NDSSName.my.fctrTStyle## 
##                                        -0.60538025 
##                                    WordCount.root2 
##                                         0.05783392 
## [1] "max lambda < lambdaOpt:"
##                                        (Intercept) 
##                                        -3.95644632 
##                       NDSSName.my.fctr#Multimedia# 
##                                        -0.07182859 
##             NDSSName.my.fctr#Opinion#RoomForDebate 
##                                        -1.30034382 
##           NDSSName.my.fctr#Opinion#ThePublicEditor 
##                                         3.29694236 
##                    NDSSName.my.fctr#U.S.#Education 
##                                        -0.53415905 
## NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness 
##                                        -0.14259759 
##         NDSSName.my.fctrBusiness#Crosswords/Games# 
##                                         3.94231812 
##               NDSSName.my.fctrBusiness#Technology# 
##                                         0.45657914 
##                     NDSSName.my.fctrForeign#World# 
##                                        -0.10021084 
##          NDSSName.my.fctrForeign#World#AsiaPacific 
##                                        -0.53077048 
##                      NDSSName.my.fctrOpEd#Opinion# 
##                                         4.01324666 
##                    NDSSName.my.fctrScience#Health# 
##                                         3.18936803 
##                    NDSSName.my.fctrStyles##Fashion 
##                                        -0.38069674 
##                       NDSSName.my.fctrStyles#U.S.# 
##                                         2.97176051 
##                           NDSSName.my.fctrTStyle## 
##                                        -0.64837249 
##                                    WordCount.root2 
##                                         0.05978318
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-26.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-27.png) 

```
##   Popular.fctr Popular.fctr.predict.Max.cor.Y.rcv.3X5.glmnet.N
## 1            N                                            3796
## 2            Y                                             177
##   Popular.fctr.predict.Max.cor.Y.rcv.3X5.glmnet.Y
## 1                                             145
## 2                                             686
##          Prediction
## Reference    N    Y
##         N 3796  145
##         Y  177  686
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   9.329725e-01   7.692476e-01   9.255302e-01   9.398832e-01   8.203580e-01 
## AccuracyPValue  McnemarPValue 
##  1.026390e-114   8.406670e-02
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-28.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-29.png) 

```
##   Popular.fctr Popular.fctr.predict.Max.cor.Y.rcv.3X5.glmnet.N
## 1            N                                            1146
## 2            Y                                              67
##   Popular.fctr.predict.Max.cor.Y.rcv.3X5.glmnet.Y
## 1                                             352
## 2                                             163
##          Prediction
## Reference    N    Y
##         N 1146  352
##         Y   67  163
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   7.575231e-01   3.107477e-01   7.365992e-01   7.775689e-01   8.668981e-01 
## AccuracyPValue  McnemarPValue 
##   1.000000e+00   9.066396e-44 
##                         id                            feats
## 1 Max.cor.Y.rcv.3X5.glmnet WordCount.root2,NDSSName.my.fctr
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1              25                      6.419                 0.269
##   max.AUCpROC.fit max.Sens.fit max.Spec.fit max.AUCROCR.fit
## 1       0.8767919     0.964476    0.7891078       0.9582555
##   opt.prob.threshold.fit max.f.score.fit max.Accuracy.fit
## 1                    0.4       0.8099174        0.9332218
##   max.AccuracyLower.fit max.AccuracyUpper.fit max.Kappa.fit
## 1             0.9255302             0.9398832     0.7686375
##   max.AUCpROC.OOB max.Sens.OOB max.Spec.OOB max.AUCROCR.OOB
## 1       0.5962443    0.9098798    0.2826087       0.8067975
##   opt.prob.threshold.OOB max.f.score.OOB max.Accuracy.OOB
## 1                    0.1       0.4375839        0.7575231
##   max.AccuracyLower.OOB max.AccuracyUpper.OOB max.Kappa.OOB
## 1             0.7365992             0.7775689     0.3107477
##   max.AccuracySD.fit max.KappaSD.fit
## 1        0.005396525      0.01835474
## [1] "fitting model: Max.cor.Y.rcv.5X1.glmnet"
## [1] "    indep_vars: WordCount.root2,NDSSName.my.fctr"
## Aggregating results
## Selecting tuning parameters
## Fitting alpha = 0.1, lambda = 0.0201 on full training set
```

```
## Warning in myfit_mdl(mdl_specs_lst = myinit_mdl_specs_lst(mdl_specs_lst =
## list(id.prefix = paste0("Max.cor.Y.rcv.", : model's bestTune found at an
## extreme of tuneGrid for parameter: alpha
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-30.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-31.png) 

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
##                                        (Intercept) 
##                                        -3.81141260 
##                       NDSSName.my.fctr#Multimedia# 
##                                        -0.68105584 
##             NDSSName.my.fctr#Opinion#RoomForDebate 
##                                        -1.92624537 
##           NDSSName.my.fctr#Opinion#ThePublicEditor 
##                                         3.40699589 
##                    NDSSName.my.fctr#U.S.#Education 
##                                        -0.98291999 
##      NDSSName.my.fctrBusiness#BusinessDay#Dealbook 
##                                        -0.22577146 
## NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness 
##                                        -0.64343834 
##         NDSSName.my.fctrBusiness#Crosswords/Games# 
##                                         3.82797332 
##               NDSSName.my.fctrBusiness#Technology# 
##                                         0.45317927 
##                      NDSSName.my.fctrCulture#Arts# 
##                                        -0.17187706 
##                     NDSSName.my.fctrForeign#World# 
##                                        -0.72035867 
##          NDSSName.my.fctrForeign#World#AsiaPacific 
##                                        -0.99018968 
##                      NDSSName.my.fctrOpEd#Opinion# 
##                                         3.81891156 
##                    NDSSName.my.fctrScience#Health# 
##                                         3.05516080 
##                    NDSSName.my.fctrStyles##Fashion 
##                                        -0.97651721 
##                       NDSSName.my.fctrStyles#U.S.# 
##                                         2.84779285 
##                           NDSSName.my.fctrTStyle## 
##                                        -0.94109645 
##                     NDSSName.my.fctrTravel#Travel# 
##                                        -0.68827560 
##                            NDSSName.my.fctrmyOther 
##                                        -0.84423735 
##                                    WordCount.root2 
##                                         0.06115867 
## [1] "max lambda < lambdaOpt:"
##                                        (Intercept) 
##                                        -3.87108412 
##                       NDSSName.my.fctr#Multimedia# 
##                                        -0.71588942 
##             NDSSName.my.fctr#Opinion#RoomForDebate 
##                                        -2.02010163 
##           NDSSName.my.fctr#Opinion#ThePublicEditor 
##                                         3.46715540 
##                    NDSSName.my.fctr#U.S.#Education 
##                                        -1.02957582 
##      NDSSName.my.fctrBusiness#BusinessDay#Dealbook 
##                                        -0.22558850 
## NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness 
##                                        -0.66798026 
##         NDSSName.my.fctrBusiness#Crosswords/Games# 
##                                         3.89132347 
##               NDSSName.my.fctrBusiness#Technology# 
##                                         0.48212450 
##                      NDSSName.my.fctrCulture#Arts# 
##                                        -0.16733777 
##                     NDSSName.my.fctrForeign#World# 
##                                        -0.75793881 
##          NDSSName.my.fctrForeign#World#AsiaPacific 
##                                        -1.03076807 
##                      NDSSName.my.fctrOpEd#Opinion# 
##                                         3.87908175 
##                    NDSSName.my.fctrScience#Health# 
##                                         3.09788786 
##                    NDSSName.my.fctrStyles##Fashion 
##                                        -1.02481879 
##                       NDSSName.my.fctrStyles#U.S.# 
##                                         2.88826078 
##                           NDSSName.my.fctrTStyle## 
##                                        -0.97585470 
##                     NDSSName.my.fctrTravel#Travel# 
##                                        -0.72668427 
##                            NDSSName.my.fctrmyOther 
##                                        -0.90347045 
##                                    WordCount.root2 
##                                         0.06289698
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-32.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-33.png) 

```
##   Popular.fctr Popular.fctr.predict.Max.cor.Y.rcv.5X1.glmnet.N
## 1            N                                            3800
## 2            Y                                             179
##   Popular.fctr.predict.Max.cor.Y.rcv.5X1.glmnet.Y
## 1                                             141
## 2                                             684
##          Prediction
## Reference    N    Y
##         N 3800  141
##         Y  179  684
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   9.333888e-01   7.700473e-01   9.259666e-01   9.402789e-01   8.203580e-01 
## AccuracyPValue  McnemarPValue 
##  1.097051e-115   3.860591e-02
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-34.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-35.png) 

```
##   Popular.fctr Popular.fctr.predict.Max.cor.Y.rcv.5X1.glmnet.N
## 1            N                                            1137
## 2            Y                                              53
##   Popular.fctr.predict.Max.cor.Y.rcv.5X1.glmnet.Y
## 1                                             361
## 2                                             177
##          Prediction
## Reference    N    Y
##         N 1137  361
##         Y   53  177
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   7.604167e-01   3.373693e-01   7.395703e-01   7.803749e-01   8.668981e-01 
## AccuracyPValue  McnemarPValue 
##   1.000000e+00   1.935747e-51 
##                         id                            feats
## 1 Max.cor.Y.rcv.5X1.glmnet WordCount.root2,NDSSName.my.fctr
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1              25                      3.223                 0.267
##   max.AUCpROC.fit max.Sens.fit max.Spec.fit max.AUCROCR.fit
## 1       0.8784031    0.9642223     0.792584       0.9607052
##   opt.prob.threshold.fit max.f.score.fit max.Accuracy.fit
## 1                    0.5       0.8104265        0.9331818
##   max.AccuracyLower.fit max.AccuracyUpper.fit max.Kappa.fit
## 1             0.9259666             0.9402789     0.7689055
##   max.AUCpROC.OOB max.Sens.OOB max.Spec.OOB max.AUCROCR.OOB
## 1       0.5962443    0.9098798    0.2826087       0.8114863
##   opt.prob.threshold.OOB max.f.score.OOB max.Accuracy.OOB
## 1                    0.1       0.4609375        0.7604167
##   max.AccuracyLower.OOB max.AccuracyUpper.OOB max.Kappa.OOB
## 1             0.7395703             0.7803749     0.3373693
##   max.AccuracySD.fit max.KappaSD.fit
## 1        0.008837283      0.03133449
## [1] "fitting model: Max.cor.Y.rcv.5X3.glmnet"
## [1] "    indep_vars: WordCount.root2,NDSSName.my.fctr"
## Aggregating results
## Selecting tuning parameters
## Fitting alpha = 0.1, lambda = 0.0201 on full training set
```

```
## Warning in myfit_mdl(mdl_specs_lst = myinit_mdl_specs_lst(mdl_specs_lst =
## list(id.prefix = paste0("Max.cor.Y.rcv.", : model's bestTune found at an
## extreme of tuneGrid for parameter: alpha
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-36.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-37.png) 

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
##                                        (Intercept) 
##                                        -3.81141260 
##                       NDSSName.my.fctr#Multimedia# 
##                                        -0.68105584 
##             NDSSName.my.fctr#Opinion#RoomForDebate 
##                                        -1.92624537 
##           NDSSName.my.fctr#Opinion#ThePublicEditor 
##                                         3.40699589 
##                    NDSSName.my.fctr#U.S.#Education 
##                                        -0.98291999 
##      NDSSName.my.fctrBusiness#BusinessDay#Dealbook 
##                                        -0.22577146 
## NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness 
##                                        -0.64343834 
##         NDSSName.my.fctrBusiness#Crosswords/Games# 
##                                         3.82797332 
##               NDSSName.my.fctrBusiness#Technology# 
##                                         0.45317927 
##                      NDSSName.my.fctrCulture#Arts# 
##                                        -0.17187706 
##                     NDSSName.my.fctrForeign#World# 
##                                        -0.72035867 
##          NDSSName.my.fctrForeign#World#AsiaPacific 
##                                        -0.99018968 
##                      NDSSName.my.fctrOpEd#Opinion# 
##                                         3.81891156 
##                    NDSSName.my.fctrScience#Health# 
##                                         3.05516080 
##                    NDSSName.my.fctrStyles##Fashion 
##                                        -0.97651721 
##                       NDSSName.my.fctrStyles#U.S.# 
##                                         2.84779285 
##                           NDSSName.my.fctrTStyle## 
##                                        -0.94109645 
##                     NDSSName.my.fctrTravel#Travel# 
##                                        -0.68827560 
##                            NDSSName.my.fctrmyOther 
##                                        -0.84423735 
##                                    WordCount.root2 
##                                         0.06115867 
## [1] "max lambda < lambdaOpt:"
##                                        (Intercept) 
##                                        -3.87108412 
##                       NDSSName.my.fctr#Multimedia# 
##                                        -0.71588942 
##             NDSSName.my.fctr#Opinion#RoomForDebate 
##                                        -2.02010163 
##           NDSSName.my.fctr#Opinion#ThePublicEditor 
##                                         3.46715540 
##                    NDSSName.my.fctr#U.S.#Education 
##                                        -1.02957582 
##      NDSSName.my.fctrBusiness#BusinessDay#Dealbook 
##                                        -0.22558850 
## NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness 
##                                        -0.66798026 
##         NDSSName.my.fctrBusiness#Crosswords/Games# 
##                                         3.89132347 
##               NDSSName.my.fctrBusiness#Technology# 
##                                         0.48212450 
##                      NDSSName.my.fctrCulture#Arts# 
##                                        -0.16733777 
##                     NDSSName.my.fctrForeign#World# 
##                                        -0.75793881 
##          NDSSName.my.fctrForeign#World#AsiaPacific 
##                                        -1.03076807 
##                      NDSSName.my.fctrOpEd#Opinion# 
##                                         3.87908175 
##                    NDSSName.my.fctrScience#Health# 
##                                         3.09788786 
##                    NDSSName.my.fctrStyles##Fashion 
##                                        -1.02481879 
##                       NDSSName.my.fctrStyles#U.S.# 
##                                         2.88826078 
##                           NDSSName.my.fctrTStyle## 
##                                        -0.97585470 
##                     NDSSName.my.fctrTravel#Travel# 
##                                        -0.72668427 
##                            NDSSName.my.fctrmyOther 
##                                        -0.90347045 
##                                    WordCount.root2 
##                                         0.06289698
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-38.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-39.png) 

```
##   Popular.fctr Popular.fctr.predict.Max.cor.Y.rcv.5X3.glmnet.N
## 1            N                                            3800
## 2            Y                                             179
##   Popular.fctr.predict.Max.cor.Y.rcv.5X3.glmnet.Y
## 1                                             141
## 2                                             684
##          Prediction
## Reference    N    Y
##         N 3800  141
##         Y  179  684
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   9.333888e-01   7.700473e-01   9.259666e-01   9.402789e-01   8.203580e-01 
## AccuracyPValue  McnemarPValue 
##  1.097051e-115   3.860591e-02
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-40.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-41.png) 

```
##   Popular.fctr Popular.fctr.predict.Max.cor.Y.rcv.5X3.glmnet.N
## 1            N                                            1137
## 2            Y                                              53
##   Popular.fctr.predict.Max.cor.Y.rcv.5X3.glmnet.Y
## 1                                             361
## 2                                             177
##          Prediction
## Reference    N    Y
##         N 1137  361
##         Y   53  177
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   7.604167e-01   3.373693e-01   7.395703e-01   7.803749e-01   8.668981e-01 
## AccuracyPValue  McnemarPValue 
##   1.000000e+00   1.935747e-51 
##                         id                            feats
## 1 Max.cor.Y.rcv.5X3.glmnet WordCount.root2,NDSSName.my.fctr
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1              25                      6.142                 0.267
##   max.AUCpROC.fit max.Sens.fit max.Spec.fit max.AUCROCR.fit
## 1       0.8784031    0.9642223     0.792584       0.9607052
##   opt.prob.threshold.fit max.f.score.fit max.Accuracy.fit
## 1                    0.5       0.8104265        0.9333905
##   max.AccuracyLower.fit max.AccuracyUpper.fit max.Kappa.fit
## 1             0.9259666             0.9402789     0.7698577
##   max.AUCpROC.OOB max.Sens.OOB max.Spec.OOB max.AUCROCR.OOB
## 1       0.5962443    0.9098798    0.2826087       0.8114863
##   opt.prob.threshold.OOB max.f.score.OOB max.Accuracy.OOB
## 1                    0.1       0.4609375        0.7604167
##   max.AccuracyLower.OOB max.AccuracyUpper.OOB max.Kappa.OOB
## 1             0.7395703             0.7803749     0.3373693
##   max.AccuracySD.fit max.KappaSD.fit
## 1        0.006138477      0.02161286
## [1] "fitting model: Max.cor.Y.rcv.5X5.glmnet"
## [1] "    indep_vars: WordCount.root2,NDSSName.my.fctr"
## Aggregating results
## Selecting tuning parameters
## Fitting alpha = 0.1, lambda = 0.0201 on full training set
```

```
## Warning in myfit_mdl(mdl_specs_lst = myinit_mdl_specs_lst(mdl_specs_lst =
## list(id.prefix = paste0("Max.cor.Y.rcv.", : model's bestTune found at an
## extreme of tuneGrid for parameter: alpha
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-42.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-43.png) 

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
##                                        (Intercept) 
##                                        -3.81141260 
##                       NDSSName.my.fctr#Multimedia# 
##                                        -0.68105584 
##             NDSSName.my.fctr#Opinion#RoomForDebate 
##                                        -1.92624537 
##           NDSSName.my.fctr#Opinion#ThePublicEditor 
##                                         3.40699589 
##                    NDSSName.my.fctr#U.S.#Education 
##                                        -0.98291999 
##      NDSSName.my.fctrBusiness#BusinessDay#Dealbook 
##                                        -0.22577146 
## NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness 
##                                        -0.64343834 
##         NDSSName.my.fctrBusiness#Crosswords/Games# 
##                                         3.82797332 
##               NDSSName.my.fctrBusiness#Technology# 
##                                         0.45317927 
##                      NDSSName.my.fctrCulture#Arts# 
##                                        -0.17187706 
##                     NDSSName.my.fctrForeign#World# 
##                                        -0.72035867 
##          NDSSName.my.fctrForeign#World#AsiaPacific 
##                                        -0.99018968 
##                      NDSSName.my.fctrOpEd#Opinion# 
##                                         3.81891156 
##                    NDSSName.my.fctrScience#Health# 
##                                         3.05516080 
##                    NDSSName.my.fctrStyles##Fashion 
##                                        -0.97651721 
##                       NDSSName.my.fctrStyles#U.S.# 
##                                         2.84779285 
##                           NDSSName.my.fctrTStyle## 
##                                        -0.94109645 
##                     NDSSName.my.fctrTravel#Travel# 
##                                        -0.68827560 
##                            NDSSName.my.fctrmyOther 
##                                        -0.84423735 
##                                    WordCount.root2 
##                                         0.06115867 
## [1] "max lambda < lambdaOpt:"
##                                        (Intercept) 
##                                        -3.87108412 
##                       NDSSName.my.fctr#Multimedia# 
##                                        -0.71588942 
##             NDSSName.my.fctr#Opinion#RoomForDebate 
##                                        -2.02010163 
##           NDSSName.my.fctr#Opinion#ThePublicEditor 
##                                         3.46715540 
##                    NDSSName.my.fctr#U.S.#Education 
##                                        -1.02957582 
##      NDSSName.my.fctrBusiness#BusinessDay#Dealbook 
##                                        -0.22558850 
## NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness 
##                                        -0.66798026 
##         NDSSName.my.fctrBusiness#Crosswords/Games# 
##                                         3.89132347 
##               NDSSName.my.fctrBusiness#Technology# 
##                                         0.48212450 
##                      NDSSName.my.fctrCulture#Arts# 
##                                        -0.16733777 
##                     NDSSName.my.fctrForeign#World# 
##                                        -0.75793881 
##          NDSSName.my.fctrForeign#World#AsiaPacific 
##                                        -1.03076807 
##                      NDSSName.my.fctrOpEd#Opinion# 
##                                         3.87908175 
##                    NDSSName.my.fctrScience#Health# 
##                                         3.09788786 
##                    NDSSName.my.fctrStyles##Fashion 
##                                        -1.02481879 
##                       NDSSName.my.fctrStyles#U.S.# 
##                                         2.88826078 
##                           NDSSName.my.fctrTStyle## 
##                                        -0.97585470 
##                     NDSSName.my.fctrTravel#Travel# 
##                                        -0.72668427 
##                            NDSSName.my.fctrmyOther 
##                                        -0.90347045 
##                                    WordCount.root2 
##                                         0.06289698
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-44.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-45.png) 

```
##   Popular.fctr Popular.fctr.predict.Max.cor.Y.rcv.5X5.glmnet.N
## 1            N                                            3800
## 2            Y                                             179
##   Popular.fctr.predict.Max.cor.Y.rcv.5X5.glmnet.Y
## 1                                             141
## 2                                             684
##          Prediction
## Reference    N    Y
##         N 3800  141
##         Y  179  684
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   9.333888e-01   7.700473e-01   9.259666e-01   9.402789e-01   8.203580e-01 
## AccuracyPValue  McnemarPValue 
##  1.097051e-115   3.860591e-02
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-46.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-47.png) 

```
##   Popular.fctr Popular.fctr.predict.Max.cor.Y.rcv.5X5.glmnet.N
## 1            N                                            1137
## 2            Y                                              53
##   Popular.fctr.predict.Max.cor.Y.rcv.5X5.glmnet.Y
## 1                                             361
## 2                                             177
##          Prediction
## Reference    N    Y
##         N 1137  361
##         Y   53  177
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   7.604167e-01   3.373693e-01   7.395703e-01   7.803749e-01   8.668981e-01 
## AccuracyPValue  McnemarPValue 
##   1.000000e+00   1.935747e-51 
##                         id                            feats
## 1 Max.cor.Y.rcv.5X5.glmnet WordCount.root2,NDSSName.my.fctr
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1              25                      8.908                 0.268
##   max.AUCpROC.fit max.Sens.fit max.Spec.fit max.AUCROCR.fit
## 1       0.8784031    0.9642223     0.792584       0.9607052
##   opt.prob.threshold.fit max.f.score.fit max.Accuracy.fit
## 1                    0.5       0.8104265        0.9331816
##   max.AccuracyLower.fit max.AccuracyUpper.fit max.Kappa.fit
## 1             0.9259666             0.9402789     0.7691429
##   max.AUCpROC.OOB max.Sens.OOB max.Spec.OOB max.AUCROCR.OOB
## 1       0.5962443    0.9098798    0.2826087       0.8114863
##   opt.prob.threshold.OOB max.f.score.OOB max.Accuracy.OOB
## 1                    0.1       0.4609375        0.7604167
##   max.AccuracyLower.OOB max.AccuracyUpper.OOB max.Kappa.OOB
## 1             0.7395703             0.7803749     0.3373693
##   max.AccuracySD.fit max.KappaSD.fit
## 1          0.0062138      0.02210061
```

```r
# Add parallel coordinates graph of glb_models_df[, glbMdlMetricsEval] to evaluate cv parameters
tmp_models_cols <- c("id", "max.nTuningRuns",
                    glbMdlMetricsEval[glbMdlMetricsEval %in% names(glb_models_df)],
                    grep("opt.", names(glb_models_df), fixed = TRUE, value = TRUE)) 
print(myplot_parcoord(obs_df = subset(glb_models_df, 
                                      grepl("Max.cor.Y.rcv.", id, fixed = TRUE), 
                                        select = -feats)[, tmp_models_cols],
                      id_var = "id"))
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-48.png) 

```r
ret_lst <- myfit_mdl(mdl_specs_lst=myinit_mdl_specs_lst(mdl_specs_lst=list(
    id.prefix="Max.cor.Y.rcv.1X1.cp.0", type=glb_model_type, trainControl.method="none",
    train.method="rpart",
    tune.df=data.frame(method="rpart", parameter="cp", min=0.0, max=0.0, by=0.1))),
                    indep_vars=max_cor_y_x_vars, rsp_var=glb_rsp_var, 
                    fit_df=glbObsFit, OOB_df=glbObsOOB)
```

```
## [1] "fitting model: Max.cor.Y.rcv.1X1.cp.0.rpart"
## [1] "    indep_vars: WordCount.root2,NDSSName.my.fctr"
```

```
## Loading required package: rpart
```

```
## Fitting cp = 0 on full training set
```

```
## Loading required package: rpart.plot
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-49.png) 

```
## Call:
## rpart(formula = .outcome ~ ., control = list(minsplit = 20, minbucket = 7, 
##     cp = 0, maxcompete = 4, maxsurrogate = 5, usesurrogate = 2, 
##     surrogatestyle = 0, maxdepth = 30, xval = 0))
##   n= 4804 
## 
##              CP nsplit rel error
## 1  0.3696407879      0 1.0000000
## 2  0.0984936269      1 0.6303592
## 3  0.0857473928      2 0.5318656
## 4  0.0567786790      3 0.4461182
## 5  0.0104287370      4 0.3893395
## 6  0.0057937428      5 0.3789108
## 7  0.0034762457      7 0.3673233
## 8  0.0023174971      8 0.3638470
## 9  0.0011587486     11 0.3568946
## 10 0.0007724990     13 0.3545771
## 11 0.0005793743     16 0.3522596
## 12 0.0004213631     24 0.3476246
## 13 0.0003862495     35 0.3429896
## 14 0.0000000000     41 0.3406721
## 
## Variable importance
##                 NDSSName.my.fctrOpEd#Opinion# 
##                                            48 
##    NDSSName.my.fctrBusiness#Crosswords/Games# 
##                                            14 
##               NDSSName.my.fctrScience#Health# 
##                                            14 
##                  NDSSName.my.fctrStyles#U.S.# 
##                                            11 
##                               WordCount.root2 
##                                             9 
##      NDSSName.my.fctr#Opinion#ThePublicEditor 
##                                             2 
## NDSSName.my.fctrBusiness#BusinessDay#Dealbook 
##                                             1 
## 
## Node number 1: 4804 observations,    complexity param=0.3696408
##   predicted class=N  expected loss=0.179642  P(node) =1
##     class counts:  3941   863
##    probabilities: 0.820 0.180 
##   left son=2 (4367 obs) right son=3 (437 obs)
##   Primary splits:
##       NDSSName.my.fctrOpEd#Opinion#              < 0.5      to the left,  improve=451.59770, (0 missing)
##       NDSSName.my.fctrBusiness#Crosswords/Games# < 0.5      to the left,  improve=112.88510, (0 missing)
##       WordCount.root2                            < 25.75849 to the left,  improve=111.17610, (0 missing)
##       NDSSName.my.fctrScience#Health#            < 0.5      to the left,  improve= 99.35206, (0 missing)
##       NDSSName.my.fctrStyles#U.S.#               < 0.5      to the left,  improve= 68.73272, (0 missing)
## 
## Node number 2: 4367 observations,    complexity param=0.09849363
##   predicted class=N  expected loss=0.1110602  P(node) =0.9090341
##     class counts:  3882   485
##    probabilities: 0.889 0.111 
##   left son=4 (4262 obs) right son=5 (105 obs)
##   Primary splits:
##       NDSSName.my.fctrBusiness#Crosswords/Games# < 0.5      to the left,  improve=135.55130, (0 missing)
##       NDSSName.my.fctrScience#Health#            < 0.5      to the left,  improve=125.07920, (0 missing)
##       WordCount.root2                            < 25.75849 to the left,  improve= 94.70710, (0 missing)
##       NDSSName.my.fctrStyles#U.S.#               < 0.5      to the left,  improve= 88.56821, (0 missing)
##       NDSSName.my.fctr#Opinion#ThePublicEditor   < 0.5      to the left,  improve= 18.74400, (0 missing)
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
##       NDSSName.my.fctrScience#Health#          < 0.5      to the left,  improve=132.96710, (0 missing)
##       NDSSName.my.fctrStyles#U.S.#             < 0.5      to the left,  improve= 94.69099, (0 missing)
##       WordCount.root2                          < 26.49528 to the left,  improve= 84.07487, (0 missing)
##       NDSSName.my.fctr#Opinion#ThePublicEditor < 0.5      to the left,  improve= 19.71762, (0 missing)
##       NDSSName.my.fctrTStyle##                 < 0.5      to the right, improve= 10.17000, (0 missing)
## 
## Node number 5: 105 observations,    complexity param=0.002317497
##   predicted class=Y  expected loss=0.0952381  P(node) =0.02185679
##     class counts:    10    95
##    probabilities: 0.095 0.905 
##   left son=10 (12 obs) right son=11 (93 obs)
##   Primary splits:
##       WordCount.root2 < 18.9043  to the left,  improve=6.455453, (0 missing)
## 
## Node number 8: 4114 observations,    complexity param=0.05677868
##   predicted class=N  expected loss=0.06781721  P(node) =0.8563697
##     class counts:  3835   279
##    probabilities: 0.932 0.068 
##   left son=16 (3987 obs) right son=17 (127 obs)
##   Primary splits:
##       NDSSName.my.fctrStyles#U.S.#             < 0.5      to the left,  improve=102.410700, (0 missing)
##       WordCount.root2                          < 25.01    to the left,  improve= 47.352210, (0 missing)
##       NDSSName.my.fctr#Opinion#ThePublicEditor < 0.5      to the left,  improve= 20.930810, (0 missing)
##       NDSSName.my.fctrTStyle##                 < 0.5      to the right, improve=  5.249425, (0 missing)
##       NDSSName.my.fctrBusiness#Technology#     < 0.5      to the left,  improve=  2.395935, (0 missing)
## 
## Node number 9: 148 observations,    complexity param=0.01042874
##   predicted class=Y  expected loss=0.25  P(node) =0.03080766
##     class counts:    37   111
##    probabilities: 0.250 0.750 
##   left son=18 (55 obs) right son=19 (93 obs)
##   Primary splits:
##       WordCount.root2 < 22.72663 to the left,  improve=19.274, (0 missing)
## 
## Node number 10: 12 observations
##   predicted class=N  expected loss=0.4166667  P(node) =0.002497918
##     class counts:     7     5
##    probabilities: 0.583 0.417 
## 
## Node number 11: 93 observations
##   predicted class=Y  expected loss=0.03225806  P(node) =0.01935887
##     class counts:     3    90
##    probabilities: 0.032 0.968 
## 
## Node number 16: 3987 observations,    complexity param=0.005793743
##   predicted class=N  expected loss=0.04790569  P(node) =0.8299334
##     class counts:  3796   191
##    probabilities: 0.952 0.048 
##   left son=32 (2982 obs) right son=33 (1005 obs)
##   Primary splits:
##       WordCount.root2                          < 25.01    to the left,  improve=29.253580, (0 missing)
##       NDSSName.my.fctr#Opinion#ThePublicEditor < 0.5      to the left,  improve=21.978920, (0 missing)
##       NDSSName.my.fctrBusiness#Technology#     < 0.5      to the left,  improve= 3.887348, (0 missing)
##       NDSSName.my.fctrTStyle##                 < 0.5      to the right, improve= 2.348653, (0 missing)
##       NDSSName.my.fctr#U.S.#Education          < 0.5      to the right, improve= 1.187739, (0 missing)
##   Surrogate splits:
##       NDSSName.my.fctr#Opinion#RoomForDebate             < 0.5      to the left,  agree=0.758, adj=0.042, (0 split)
##       NDSSName.my.fctrForeign#World#AsiaPacific          < 0.5      to the left,  agree=0.752, adj=0.016, (0 split)
##       NDSSName.my.fctr#Opinion#ThePublicEditor           < 0.5      to the left,  agree=0.750, adj=0.008, (0 split)
##       NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness < 0.5      to the left,  agree=0.748, adj=0.002, (0 split)
## 
## Node number 17: 127 observations,    complexity param=0.003476246
##   predicted class=Y  expected loss=0.3070866  P(node) =0.0264363
##     class counts:    39    88
##    probabilities: 0.307 0.693 
##   left son=34 (13 obs) right son=35 (114 obs)
##   Primary splits:
##       WordCount.root2 < 15.32846 to the left,  improve=2.753047, (0 missing)
## 
## Node number 18: 55 observations,    complexity param=0.002317497
##   predicted class=N  expected loss=0.4181818  P(node) =0.01144879
##     class counts:    32    23
##    probabilities: 0.582 0.418 
##   left son=36 (9 obs) right son=37 (46 obs)
##   Primary splits:
##       WordCount.root2 < 19.93708 to the right, improve=0.8264383, (0 missing)
## 
## Node number 19: 93 observations
##   predicted class=Y  expected loss=0.05376344  P(node) =0.01935887
##     class counts:     5    88
##    probabilities: 0.054 0.946 
## 
## Node number 32: 2982 observations
##   predicted class=N  expected loss=0.01274313  P(node) =0.6207327
##     class counts:  2944    38
##    probabilities: 0.987 0.013 
## 
## Node number 33: 1005 observations,    complexity param=0.005793743
##   predicted class=N  expected loss=0.1522388  P(node) =0.2092007
##     class counts:   852   153
##    probabilities: 0.848 0.152 
##   left son=66 (993 obs) right son=67 (12 obs)
##   Primary splits:
##       NDSSName.my.fctr#Opinion#ThePublicEditor  < 0.5      to the left,  improve=14.193880, (0 missing)
##       NDSSName.my.fctrCulture#Arts#             < 0.5      to the left,  improve= 3.669601, (0 missing)
##       NDSSName.my.fctrForeign#World#AsiaPacific < 0.5      to the right, improve= 3.556158, (0 missing)
##       WordCount.root2                           < 34.19795 to the left,  improve= 2.582851, (0 missing)
##       NDSSName.my.fctr#Opinion#RoomForDebate    < 0.5      to the right, improve= 2.031748, (0 missing)
## 
## Node number 34: 13 observations
##   predicted class=N  expected loss=0.3846154  P(node) =0.002706078
##     class counts:     8     5
##    probabilities: 0.615 0.385 
## 
## Node number 35: 114 observations,    complexity param=0.000772499
##   predicted class=Y  expected loss=0.2719298  P(node) =0.02373022
##     class counts:    31    83
##    probabilities: 0.272 0.728 
##   left son=70 (79 obs) right son=71 (35 obs)
##   Primary splits:
##       WordCount.root2 < 29.21444 to the left,  improve=1.020279, (0 missing)
## 
## Node number 36: 9 observations
##   predicted class=N  expected loss=0.2222222  P(node) =0.001873439
##     class counts:     7     2
##    probabilities: 0.778 0.222 
## 
## Node number 37: 46 observations,    complexity param=0.002317497
##   predicted class=N  expected loss=0.4565217  P(node) =0.009575354
##     class counts:    25    21
##    probabilities: 0.543 0.457 
##   left son=74 (36 obs) right son=75 (10 obs)
##   Primary splits:
##       WordCount.root2 < 17.01454 to the left,  improve=1.514976, (0 missing)
## 
## Node number 66: 993 observations,    complexity param=0.0005793743
##   predicted class=N  expected loss=0.143001  P(node) =0.2067027
##     class counts:   851   142
##    probabilities: 0.857 0.143 
##   left son=132 (930 obs) right son=133 (63 obs)
##   Primary splits:
##       NDSSName.my.fctrCulture#Arts#             < 0.5      to the left,  improve=4.094729, (0 missing)
##       NDSSName.my.fctrForeign#World#AsiaPacific < 0.5      to the right, improve=3.106316, (0 missing)
##       WordCount.root2                           < 29.5127  to the left,  improve=2.722793, (0 missing)
##       NDSSName.my.fctrBusiness#Technology#      < 0.5      to the left,  improve=1.962300, (0 missing)
##       NDSSName.my.fctr#Opinion#RoomForDebate    < 0.5      to the right, improve=1.793603, (0 missing)
## 
## Node number 67: 12 observations
##   predicted class=Y  expected loss=0.08333333  P(node) =0.002497918
##     class counts:     1    11
##    probabilities: 0.083 0.917 
## 
## Node number 70: 79 observations,    complexity param=0.000772499
##   predicted class=Y  expected loss=0.3164557  P(node) =0.01644463
##     class counts:    25    54
##    probabilities: 0.316 0.684 
##   left son=140 (25 obs) right son=141 (54 obs)
##   Primary splits:
##       WordCount.root2 < 27.36786 to the right, improve=0.5105485, (0 missing)
## 
## Node number 71: 35 observations
##   predicted class=Y  expected loss=0.1714286  P(node) =0.007285595
##     class counts:     6    29
##    probabilities: 0.171 0.829 
## 
## Node number 74: 36 observations,    complexity param=0.001158749
##   predicted class=N  expected loss=0.3888889  P(node) =0.007493755
##     class counts:    22    14
##    probabilities: 0.611 0.389 
##   left son=148 (8 obs) right son=149 (28 obs)
##   Primary splits:
##       WordCount.root2 < 15.74773 to the right, improve=0.3968254, (0 missing)
## 
## Node number 75: 10 observations
##   predicted class=Y  expected loss=0.3  P(node) =0.002081599
##     class counts:     3     7
##    probabilities: 0.300 0.700 
## 
## Node number 132: 930 observations,    complexity param=0.0005793743
##   predicted class=N  expected loss=0.1311828  P(node) =0.1935887
##     class counts:   808   122
##    probabilities: 0.869 0.131 
##   left son=264 (627 obs) right son=265 (303 obs)
##   Primary splits:
##       WordCount.root2                               < 33.97057 to the left,  improve=2.913816, (0 missing)
##       NDSSName.my.fctrForeign#World#AsiaPacific     < 0.5      to the right, improve=2.586923, (0 missing)
##       NDSSName.my.fctrBusiness#Technology#          < 0.5      to the left,  improve=2.402029, (0 missing)
##       NDSSName.my.fctr#Opinion#RoomForDebate        < 0.5      to the right, improve=1.513920, (0 missing)
##       NDSSName.my.fctrBusiness#BusinessDay#Dealbook < 0.5      to the left,  improve=1.276783, (0 missing)
##   Surrogate splits:
##       NDSSName.my.fctr#Opinion#RoomForDebate < 0.5      to the left,  agree=0.719, adj=0.139, (0 split)
## 
## Node number 133: 63 observations,    complexity param=0.0003862495
##   predicted class=N  expected loss=0.3174603  P(node) =0.01311407
##     class counts:    43    20
##    probabilities: 0.683 0.317 
##   left son=266 (14 obs) right son=267 (49 obs)
##   Primary splits:
##       WordCount.root2 < 26.99984 to the left,  improve=0.38322, (0 missing)
## 
## Node number 140: 25 observations,    complexity param=0.000772499
##   predicted class=Y  expected loss=0.4  P(node) =0.005203997
##     class counts:    10    15
##    probabilities: 0.400 0.600 
##   left son=280 (8 obs) right son=281 (17 obs)
##   Primary splits:
##       WordCount.root2 < 28.02674 to the left,  improve=1.191176, (0 missing)
## 
## Node number 141: 54 observations,    complexity param=0.0003862495
##   predicted class=Y  expected loss=0.2777778  P(node) =0.01124063
##     class counts:    15    39
##    probabilities: 0.278 0.722 
##   left son=282 (45 obs) right son=283 (9 obs)
##   Primary splits:
##       WordCount.root2 < 26.55173 to the left,  improve=0.6, (0 missing)
## 
## Node number 148: 8 observations
##   predicted class=N  expected loss=0.25  P(node) =0.001665279
##     class counts:     6     2
##    probabilities: 0.750 0.250 
## 
## Node number 149: 28 observations,    complexity param=0.001158749
##   predicted class=N  expected loss=0.4285714  P(node) =0.005828476
##     class counts:    16    12
##    probabilities: 0.571 0.429 
##   left son=298 (20 obs) right son=299 (8 obs)
##   Primary splits:
##       WordCount.root2 < 15.06648 to the left,  improve=0.8642857, (0 missing)
## 
## Node number 264: 627 observations,    complexity param=0.0004213631
##   predicted class=N  expected loss=0.1036683  P(node) =0.1305162
##     class counts:   562    65
##    probabilities: 0.896 0.104 
##   left son=528 (561 obs) right son=529 (66 obs)
##   Primary splits:
##       NDSSName.my.fctrBusiness#Technology#      < 0.5      to the left,  improve=2.8404170, (0 missing)
##       NDSSName.my.fctrForeign#World#AsiaPacific < 0.5      to the right, improve=1.0796950, (0 missing)
##       NDSSName.my.fctrTStyle##                  < 0.5      to the right, improve=1.0670160, (0 missing)
##       WordCount.root2                           < 29.5127  to the left,  improve=0.8966879, (0 missing)
##       NDSSName.my.fctr#Multimedia#              < 0.5      to the right, improve=0.4399337, (0 missing)
## 
## Node number 265: 303 observations,    complexity param=0.0005793743
##   predicted class=N  expected loss=0.1881188  P(node) =0.06307244
##     class counts:   246    57
##    probabilities: 0.812 0.188 
##   left son=530 (222 obs) right son=531 (81 obs)
##   Primary splits:
##       NDSSName.my.fctrBusiness#BusinessDay#Dealbook < 0.5      to the left,  improve=5.4890570, (0 missing)
##       WordCount.root2                               < 38.17067 to the right, improve=5.0156320, (0 missing)
##       NDSSName.my.fctr#Opinion#RoomForDebate        < 0.5      to the right, improve=3.4510070, (0 missing)
##       NDSSName.my.fctrForeign#World#AsiaPacific     < 0.5      to the right, improve=1.5155860, (0 missing)
##       NDSSName.my.fctr#U.S.#Education               < 0.5      to the right, improve=0.8078801, (0 missing)
##   Surrogate splits:
##       WordCount.root2 < 34.08078 to the right, agree=0.739, adj=0.025, (0 split)
## 
## Node number 266: 14 observations
##   predicted class=N  expected loss=0.2142857  P(node) =0.002914238
##     class counts:    11     3
##    probabilities: 0.786 0.214 
## 
## Node number 267: 49 observations,    complexity param=0.0003862495
##   predicted class=N  expected loss=0.3469388  P(node) =0.01019983
##     class counts:    32    17
##    probabilities: 0.653 0.347 
##   left son=534 (10 obs) right son=535 (39 obs)
##   Primary splits:
##       WordCount.root2 < 41.56249 to the right, improve=0.5425432, (0 missing)
## 
## Node number 280: 8 observations
##   predicted class=N  expected loss=0.375  P(node) =0.001665279
##     class counts:     5     3
##    probabilities: 0.625 0.375 
## 
## Node number 281: 17 observations
##   predicted class=Y  expected loss=0.2941176  P(node) =0.003538718
##     class counts:     5    12
##    probabilities: 0.294 0.706 
## 
## Node number 282: 45 observations,    complexity param=0.0003862495
##   predicted class=Y  expected loss=0.3111111  P(node) =0.009367194
##     class counts:    14    31
##    probabilities: 0.311 0.689 
##   left son=564 (23 obs) right son=565 (22 obs)
##   Primary splits:
##       WordCount.root2 < 21.70252 to the right, improve=0.6050944, (0 missing)
## 
## Node number 283: 9 observations
##   predicted class=Y  expected loss=0.1111111  P(node) =0.001873439
##     class counts:     1     8
##    probabilities: 0.111 0.889 
## 
## Node number 298: 20 observations
##   predicted class=N  expected loss=0.35  P(node) =0.004163197
##     class counts:    13     7
##    probabilities: 0.650 0.350 
## 
## Node number 299: 8 observations
##   predicted class=Y  expected loss=0.375  P(node) =0.001665279
##     class counts:     3     5
##    probabilities: 0.375 0.625 
## 
## Node number 528: 561 observations,    complexity param=0.0004213631
##   predicted class=N  expected loss=0.08734403  P(node) =0.1167777
##     class counts:   512    49
##    probabilities: 0.913 0.087 
##   left son=1056 (281 obs) right son=1057 (280 obs)
##   Primary splits:
##       WordCount.root2                           < 29.33428 to the left,  improve=1.5853030, (0 missing)
##       NDSSName.my.fctrTStyle##                  < 0.5      to the right, improve=0.7645570, (0 missing)
##       NDSSName.my.fctrForeign#World#AsiaPacific < 0.5      to the right, improve=0.7250433, (0 missing)
##       NDSSName.my.fctrStyles##Fashion           < 0.5      to the right, improve=0.3000638, (0 missing)
##       NDSSName.my.fctr#Multimedia#              < 0.5      to the right, improve=0.2729836, (0 missing)
##   Surrogate splits:
##       NDSSName.my.fctrMetro#N.Y./Region#                 < 0.5      to the left,  agree=0.560, adj=0.118, (0 split)
##       NDSSName.my.fctrForeign#World#AsiaPacific          < 0.5      to the right, agree=0.533, adj=0.064, (0 split)
##       NDSSName.my.fctr#Multimedia#                       < 0.5      to the right, agree=0.524, adj=0.046, (0 split)
##       NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness < 0.5      to the left,  agree=0.515, adj=0.029, (0 split)
##       NDSSName.my.fctrStyles##Fashion                    < 0.5      to the right, agree=0.512, adj=0.021, (0 split)
## 
## Node number 529: 66 observations,    complexity param=0.0004213631
##   predicted class=N  expected loss=0.2424242  P(node) =0.01373855
##     class counts:    50    16
##    probabilities: 0.758 0.242 
##   left son=1058 (38 obs) right son=1059 (28 obs)
##   Primary splits:
##       WordCount.root2 < 27.86575 to the left,  improve=0.6070859, (0 missing)
## 
## Node number 530: 222 observations
##   predicted class=N  expected loss=0.1306306  P(node) =0.04621149
##     class counts:   193    29
##    probabilities: 0.869 0.131 
## 
## Node number 531: 81 observations,    complexity param=0.0005793743
##   predicted class=N  expected loss=0.345679  P(node) =0.01686095
##     class counts:    53    28
##    probabilities: 0.654 0.346 
##   left son=1062 (15 obs) right son=1063 (66 obs)
##   Primary splits:
##       WordCount.root2 < 41.59766 to the right, improve=2.866218, (0 missing)
## 
## Node number 534: 10 observations
##   predicted class=N  expected loss=0.2  P(node) =0.002081599
##     class counts:     8     2
##    probabilities: 0.800 0.200 
## 
## Node number 535: 39 observations,    complexity param=0.0003862495
##   predicted class=N  expected loss=0.3846154  P(node) =0.008118235
##     class counts:    24    15
##    probabilities: 0.615 0.385 
##   left son=1070 (32 obs) right son=1071 (7 obs)
##   Primary splits:
##       WordCount.root2 < 34.23387 to the left,  improve=0.595467, (0 missing)
## 
## Node number 564: 23 observations,    complexity param=0.0003862495
##   predicted class=Y  expected loss=0.3913043  P(node) =0.004787677
##     class counts:     9    14
##    probabilities: 0.391 0.609 
##   left son=1128 (7 obs) right son=1129 (16 obs)
##   Primary splits:
##       WordCount.root2 < 23.6326  to the left,  improve=0.6529503, (0 missing)
## 
## Node number 565: 22 observations
##   predicted class=Y  expected loss=0.2272727  P(node) =0.004579517
##     class counts:     5    17
##    probabilities: 0.227 0.773 
## 
## Node number 1056: 281 observations
##   predicted class=N  expected loss=0.04982206  P(node) =0.05849292
##     class counts:   267    14
##    probabilities: 0.950 0.050 
## 
## Node number 1057: 280 observations,    complexity param=0.0004213631
##   predicted class=N  expected loss=0.125  P(node) =0.05828476
##     class counts:   245    35
##    probabilities: 0.875 0.125 
##   left son=2114 (71 obs) right son=2115 (209 obs)
##   Primary splits:
##       WordCount.root2                           < 32.57299 to the right, improve=0.8968765, (0 missing)
##       NDSSName.my.fctrTStyle##                  < 0.5      to the right, improve=0.7830739, (0 missing)
##       NDSSName.my.fctrMetro#N.Y./Region#        < 0.5      to the right, improve=0.3683673, (0 missing)
##       NDSSName.my.fctr#Multimedia#              < 0.5      to the right, improve=0.3578067, (0 missing)
##       NDSSName.my.fctrForeign#World#AsiaPacific < 0.5      to the right, improve=0.3021494, (0 missing)
## 
## Node number 1058: 38 observations
##   predicted class=N  expected loss=0.1842105  P(node) =0.007910075
##     class counts:    31     7
##    probabilities: 0.816 0.184 
## 
## Node number 1059: 28 observations,    complexity param=0.0004213631
##   predicted class=N  expected loss=0.3214286  P(node) =0.005828476
##     class counts:    19     9
##    probabilities: 0.679 0.321 
##   left son=2118 (19 obs) right son=2119 (9 obs)
##   Primary splits:
##       WordCount.root2 < 28.6269  to the right, improve=1.454052, (0 missing)
## 
## Node number 1062: 15 observations
##   predicted class=N  expected loss=0.06666667  P(node) =0.003122398
##     class counts:    14     1
##    probabilities: 0.933 0.067 
## 
## Node number 1063: 66 observations,    complexity param=0.0005793743
##   predicted class=N  expected loss=0.4090909  P(node) =0.01373855
##     class counts:    39    27
##    probabilities: 0.591 0.409 
##   left son=2126 (25 obs) right son=2127 (41 obs)
##   Primary splits:
##       WordCount.root2 < 35.6581  to the left,  improve=1.341286, (0 missing)
## 
## Node number 1070: 32 observations
##   predicted class=N  expected loss=0.34375  P(node) =0.006661116
##     class counts:    21    11
##    probabilities: 0.656 0.344 
## 
## Node number 1071: 7 observations
##   predicted class=Y  expected loss=0.4285714  P(node) =0.001457119
##     class counts:     3     4
##    probabilities: 0.429 0.571 
## 
## Node number 1128: 7 observations
##   predicted class=N  expected loss=0.4285714  P(node) =0.001457119
##     class counts:     4     3
##    probabilities: 0.571 0.429 
## 
## Node number 1129: 16 observations
##   predicted class=Y  expected loss=0.3125  P(node) =0.003330558
##     class counts:     5    11
##    probabilities: 0.312 0.688 
## 
## Node number 2114: 71 observations
##   predicted class=N  expected loss=0.05633803  P(node) =0.01477935
##     class counts:    67     4
##    probabilities: 0.944 0.056 
## 
## Node number 2115: 209 observations,    complexity param=0.0004213631
##   predicted class=N  expected loss=0.1483254  P(node) =0.04350541
##     class counts:   178    31
##    probabilities: 0.852 0.148 
##   left son=4230 (12 obs) right son=4231 (197 obs)
##   Primary splits:
##       NDSSName.my.fctrTStyle##                  < 0.5      to the right, improve=0.5601729, (0 missing)
##       NDSSName.my.fctr#Multimedia#              < 0.5      to the right, improve=0.5108985, (0 missing)
##       WordCount.root2                           < 30.09153 to the right, improve=0.4980706, (0 missing)
##       NDSSName.my.fctrMetro#N.Y./Region#        < 0.5      to the right, improve=0.4241343, (0 missing)
##       NDSSName.my.fctrForeign#World#AsiaPacific < 0.5      to the right, improve=0.3390226, (0 missing)
## 
## Node number 2118: 19 observations
##   predicted class=N  expected loss=0.2105263  P(node) =0.003955037
##     class counts:    15     4
##    probabilities: 0.789 0.211 
## 
## Node number 2119: 9 observations
##   predicted class=Y  expected loss=0.4444444  P(node) =0.001873439
##     class counts:     4     5
##    probabilities: 0.444 0.556 
## 
## Node number 2126: 25 observations
##   predicted class=N  expected loss=0.28  P(node) =0.005203997
##     class counts:    18     7
##    probabilities: 0.720 0.280 
## 
## Node number 2127: 41 observations,    complexity param=0.0005793743
##   predicted class=N  expected loss=0.4878049  P(node) =0.008534555
##     class counts:    21    20
##    probabilities: 0.512 0.488 
##   left son=4254 (30 obs) right son=4255 (11 obs)
##   Primary splits:
##       WordCount.root2 < 36.31791 to the right, improve=0.6635625, (0 missing)
## 
## Node number 4230: 12 observations
##   predicted class=N  expected loss=0  P(node) =0.002497918
##     class counts:    12     0
##    probabilities: 1.000 0.000 
## 
## Node number 4231: 197 observations,    complexity param=0.0004213631
##   predicted class=N  expected loss=0.1573604  P(node) =0.04100749
##     class counts:   166    31
##    probabilities: 0.843 0.157 
##   left son=8462 (11 obs) right son=8463 (186 obs)
##   Primary splits:
##       NDSSName.my.fctr#Multimedia#                       < 0.5      to the right, improve=0.5769882, (0 missing)
##       NDSSName.my.fctrMetro#N.Y./Region#                 < 0.5      to the right, improve=0.5314217, (0 missing)
##       WordCount.root2                                    < 30.09153 to the right, improve=0.4682049, (0 missing)
##       NDSSName.my.fctrForeign#World#AsiaPacific          < 0.5      to the right, improve=0.4106319, (0 missing)
##       NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness < 0.5      to the right, improve=0.1814254, (0 missing)
## 
## Node number 4254: 30 observations,    complexity param=0.0005793743
##   predicted class=N  expected loss=0.4333333  P(node) =0.006244796
##     class counts:    17    13
##    probabilities: 0.567 0.433 
##   left son=8508 (7 obs) right son=8509 (23 obs)
##   Primary splits:
##       WordCount.root2 < 37.14159 to the left,  improve=0.3979296, (0 missing)
## 
## Node number 4255: 11 observations
##   predicted class=Y  expected loss=0.3636364  P(node) =0.002289759
##     class counts:     4     7
##    probabilities: 0.364 0.636 
## 
## Node number 8462: 11 observations
##   predicted class=N  expected loss=0  P(node) =0.002289759
##     class counts:    11     0
##    probabilities: 1.000 0.000 
## 
## Node number 8463: 186 observations,    complexity param=0.0004213631
##   predicted class=N  expected loss=0.1666667  P(node) =0.03871774
##     class counts:   155    31
##    probabilities: 0.833 0.167 
##   left son=16926 (29 obs) right son=16927 (157 obs)
##   Primary splits:
##       NDSSName.my.fctrMetro#N.Y./Region#                 < 0.5      to the right, improve=0.6559045, (0 missing)
##       NDSSName.my.fctrForeign#World#AsiaPacific          < 0.5      to the right, improve=0.4920635, (0 missing)
##       WordCount.root2                                    < 30.09153 to the right, improve=0.3890196, (0 missing)
##       NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness < 0.5      to the right, improve=0.2415584, (0 missing)
##       NDSSName.my.fctrBusiness#BusinessDay#Dealbook      < 0.5      to the left,  improve=0.0126479, (0 missing)
## 
## Node number 8508: 7 observations
##   predicted class=N  expected loss=0.2857143  P(node) =0.001457119
##     class counts:     5     2
##    probabilities: 0.714 0.286 
## 
## Node number 8509: 23 observations,    complexity param=0.0005793743
##   predicted class=N  expected loss=0.4782609  P(node) =0.004787677
##     class counts:    12    11
##    probabilities: 0.522 0.478 
##   left son=17018 (8 obs) right son=17019 (15 obs)
##   Primary splits:
##       WordCount.root2 < 38.57459 to the right, improve=0.2615942, (0 missing)
## 
## Node number 16926: 29 observations
##   predicted class=N  expected loss=0.06896552  P(node) =0.006036636
##     class counts:    27     2
##    probabilities: 0.931 0.069 
## 
## Node number 16927: 157 observations,    complexity param=0.0004213631
##   predicted class=N  expected loss=0.1847134  P(node) =0.0326811
##     class counts:   128    29
##    probabilities: 0.815 0.185 
##   left son=33854 (18 obs) right son=33855 (139 obs)
##   Primary splits:
##       NDSSName.my.fctrForeign#World#AsiaPacific          < 0.5      to the right, improve=0.67831090, (0 missing)
##       WordCount.root2                                    < 32.38827 to the left,  improve=0.61044970, (0 missing)
##       NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness < 0.5      to the right, improve=0.38816480, (0 missing)
##       NDSSName.my.fctrBusiness#BusinessDay#Dealbook      < 0.5      to the right, improve=0.01539613, (0 missing)
##   Surrogate splits:
##       WordCount.root2 < 32.51922 to the right, agree=0.892, adj=0.056, (0 split)
## 
## Node number 17018: 8 observations
##   predicted class=N  expected loss=0.375  P(node) =0.001665279
##     class counts:     5     3
##    probabilities: 0.625 0.375 
## 
## Node number 17019: 15 observations
##   predicted class=Y  expected loss=0.4666667  P(node) =0.003122398
##     class counts:     7     8
##    probabilities: 0.467 0.533 
## 
## Node number 33854: 18 observations
##   predicted class=N  expected loss=0.05555556  P(node) =0.003746878
##     class counts:    17     1
##    probabilities: 0.944 0.056 
## 
## Node number 33855: 139 observations,    complexity param=0.0004213631
##   predicted class=N  expected loss=0.2014388  P(node) =0.02893422
##     class counts:   111    28
##    probabilities: 0.799 0.201 
##   left son=67710 (102 obs) right son=67711 (37 obs)
##   Primary splits:
##       WordCount.root2                                    < 30.09153 to the right, improve=0.9266317, (0 missing)
##       NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness < 0.5      to the right, improve=0.5580040, (0 missing)
##       NDSSName.my.fctrBusiness#BusinessDay#Dealbook      < 0.5      to the right, improve=0.1306354, (0 missing)
## 
## Node number 67710: 102 observations
##   predicted class=N  expected loss=0.1666667  P(node) =0.02123231
##     class counts:    85    17
##    probabilities: 0.833 0.167 
## 
## Node number 67711: 37 observations,    complexity param=0.0004213631
##   predicted class=N  expected loss=0.2972973  P(node) =0.007701915
##     class counts:    26    11
##    probabilities: 0.703 0.297 
##   left son=135422 (30 obs) right son=135423 (7 obs)
##   Primary splits:
##       WordCount.root2                               < 29.92488 to the left,  improve=3.00231700, (0 missing)
##       NDSSName.my.fctrBusiness#BusinessDay#Dealbook < 0.5      to the left,  improve=0.01303089, (0 missing)
## 
## Node number 135422: 30 observations
##   predicted class=N  expected loss=0.2  P(node) =0.006244796
##     class counts:    24     6
##    probabilities: 0.800 0.200 
## 
## Node number 135423: 7 observations
##   predicted class=Y  expected loss=0.2857143  P(node) =0.001457119
##     class counts:     2     5
##    probabilities: 0.286 0.714 
## 
## n= 4804 
## 
## node), split, n, loss, yval, (yprob)
##       * denotes terminal node
## 
##      1) root 4804 863 N (0.82035803 0.17964197)  
##        2) NDSSName.my.fctrOpEd#Opinion#< 0.5 4367 485 N (0.88893978 0.11106022)  
##          4) NDSSName.my.fctrBusiness#Crosswords/Games#< 0.5 4262 390 N (0.90849366 0.09150634)  
##            8) NDSSName.my.fctrScience#Health#< 0.5 4114 279 N (0.93218279 0.06781721)  
##             16) NDSSName.my.fctrStyles#U.S.#< 0.5 3987 191 N (0.95209431 0.04790569)  
##               32) WordCount.root2< 25.01 2982  38 N (0.98725687 0.01274313) *
##               33) WordCount.root2>=25.01 1005 153 N (0.84776119 0.15223881)  
##                 66) NDSSName.my.fctr#Opinion#ThePublicEditor< 0.5 993 142 N (0.85699899 0.14300101)  
##                  132) NDSSName.my.fctrCulture#Arts#< 0.5 930 122 N (0.86881720 0.13118280)  
##                    264) WordCount.root2< 33.97057 627  65 N (0.89633174 0.10366826)  
##                      528) NDSSName.my.fctrBusiness#Technology#< 0.5 561  49 N (0.91265597 0.08734403)  
##                       1056) WordCount.root2< 29.33428 281  14 N (0.95017794 0.04982206) *
##                       1057) WordCount.root2>=29.33428 280  35 N (0.87500000 0.12500000)  
##                         2114) WordCount.root2>=32.57299 71   4 N (0.94366197 0.05633803) *
##                         2115) WordCount.root2< 32.57299 209  31 N (0.85167464 0.14832536)  
##                           4230) NDSSName.my.fctrTStyle##>=0.5 12   0 N (1.00000000 0.00000000) *
##                           4231) NDSSName.my.fctrTStyle##< 0.5 197  31 N (0.84263959 0.15736041)  
##                             8462) NDSSName.my.fctr#Multimedia#>=0.5 11   0 N (1.00000000 0.00000000) *
##                             8463) NDSSName.my.fctr#Multimedia#< 0.5 186  31 N (0.83333333 0.16666667)  
##                              16926) NDSSName.my.fctrMetro#N.Y./Region#>=0.5 29   2 N (0.93103448 0.06896552) *
##                              16927) NDSSName.my.fctrMetro#N.Y./Region#< 0.5 157  29 N (0.81528662 0.18471338)  
##                                33854) NDSSName.my.fctrForeign#World#AsiaPacific>=0.5 18   1 N (0.94444444 0.05555556) *
##                                33855) NDSSName.my.fctrForeign#World#AsiaPacific< 0.5 139  28 N (0.79856115 0.20143885)  
##                                  67710) WordCount.root2>=30.09153 102  17 N (0.83333333 0.16666667) *
##                                  67711) WordCount.root2< 30.09153 37  11 N (0.70270270 0.29729730)  
##                                   135422) WordCount.root2< 29.92488 30   6 N (0.80000000 0.20000000) *
##                                   135423) WordCount.root2>=29.92488 7   2 Y (0.28571429 0.71428571) *
##                      529) NDSSName.my.fctrBusiness#Technology#>=0.5 66  16 N (0.75757576 0.24242424)  
##                       1058) WordCount.root2< 27.86575 38   7 N (0.81578947 0.18421053) *
##                       1059) WordCount.root2>=27.86575 28   9 N (0.67857143 0.32142857)  
##                         2118) WordCount.root2>=28.6269 19   4 N (0.78947368 0.21052632) *
##                         2119) WordCount.root2< 28.6269 9   4 Y (0.44444444 0.55555556) *
##                    265) WordCount.root2>=33.97057 303  57 N (0.81188119 0.18811881)  
##                      530) NDSSName.my.fctrBusiness#BusinessDay#Dealbook< 0.5 222  29 N (0.86936937 0.13063063) *
##                      531) NDSSName.my.fctrBusiness#BusinessDay#Dealbook>=0.5 81  28 N (0.65432099 0.34567901)  
##                       1062) WordCount.root2>=41.59766 15   1 N (0.93333333 0.06666667) *
##                       1063) WordCount.root2< 41.59766 66  27 N (0.59090909 0.40909091)  
##                         2126) WordCount.root2< 35.6581 25   7 N (0.72000000 0.28000000) *
##                         2127) WordCount.root2>=35.6581 41  20 N (0.51219512 0.48780488)  
##                           4254) WordCount.root2>=36.31791 30  13 N (0.56666667 0.43333333)  
##                             8508) WordCount.root2< 37.14159 7   2 N (0.71428571 0.28571429) *
##                             8509) WordCount.root2>=37.14159 23  11 N (0.52173913 0.47826087)  
##                              17018) WordCount.root2>=38.57459 8   3 N (0.62500000 0.37500000) *
##                              17019) WordCount.root2< 38.57459 15   7 Y (0.46666667 0.53333333) *
##                           4255) WordCount.root2< 36.31791 11   4 Y (0.36363636 0.63636364) *
##                  133) NDSSName.my.fctrCulture#Arts#>=0.5 63  20 N (0.68253968 0.31746032)  
##                    266) WordCount.root2< 26.99984 14   3 N (0.78571429 0.21428571) *
##                    267) WordCount.root2>=26.99984 49  17 N (0.65306122 0.34693878)  
##                      534) WordCount.root2>=41.56249 10   2 N (0.80000000 0.20000000) *
##                      535) WordCount.root2< 41.56249 39  15 N (0.61538462 0.38461538)  
##                       1070) WordCount.root2< 34.23387 32  11 N (0.65625000 0.34375000) *
##                       1071) WordCount.root2>=34.23387 7   3 Y (0.42857143 0.57142857) *
##                 67) NDSSName.my.fctr#Opinion#ThePublicEditor>=0.5 12   1 Y (0.08333333 0.91666667) *
##             17) NDSSName.my.fctrStyles#U.S.#>=0.5 127  39 Y (0.30708661 0.69291339)  
##               34) WordCount.root2< 15.32846 13   5 N (0.61538462 0.38461538) *
##               35) WordCount.root2>=15.32846 114  31 Y (0.27192982 0.72807018)  
##                 70) WordCount.root2< 29.21444 79  25 Y (0.31645570 0.68354430)  
##                  140) WordCount.root2>=27.36786 25  10 Y (0.40000000 0.60000000)  
##                    280) WordCount.root2< 28.02674 8   3 N (0.62500000 0.37500000) *
##                    281) WordCount.root2>=28.02674 17   5 Y (0.29411765 0.70588235) *
##                  141) WordCount.root2< 27.36786 54  15 Y (0.27777778 0.72222222)  
##                    282) WordCount.root2< 26.55173 45  14 Y (0.31111111 0.68888889)  
##                      564) WordCount.root2>=21.70252 23   9 Y (0.39130435 0.60869565)  
##                       1128) WordCount.root2< 23.6326 7   3 N (0.57142857 0.42857143) *
##                       1129) WordCount.root2>=23.6326 16   5 Y (0.31250000 0.68750000) *
##                      565) WordCount.root2< 21.70252 22   5 Y (0.22727273 0.77272727) *
##                    283) WordCount.root2>=26.55173 9   1 Y (0.11111111 0.88888889) *
##                 71) WordCount.root2>=29.21444 35   6 Y (0.17142857 0.82857143) *
##            9) NDSSName.my.fctrScience#Health#>=0.5 148  37 Y (0.25000000 0.75000000)  
##             18) WordCount.root2< 22.72663 55  23 N (0.58181818 0.41818182)  
##               36) WordCount.root2>=19.93708 9   2 N (0.77777778 0.22222222) *
##               37) WordCount.root2< 19.93708 46  21 N (0.54347826 0.45652174)  
##                 74) WordCount.root2< 17.01454 36  14 N (0.61111111 0.38888889)  
##                  148) WordCount.root2>=15.74773 8   2 N (0.75000000 0.25000000) *
##                  149) WordCount.root2< 15.74773 28  12 N (0.57142857 0.42857143)  
##                    298) WordCount.root2< 15.06648 20   7 N (0.65000000 0.35000000) *
##                    299) WordCount.root2>=15.06648 8   3 Y (0.37500000 0.62500000) *
##                 75) WordCount.root2>=17.01454 10   3 Y (0.30000000 0.70000000) *
##             19) WordCount.root2>=22.72663 93   5 Y (0.05376344 0.94623656) *
##          5) NDSSName.my.fctrBusiness#Crosswords/Games#>=0.5 105  10 Y (0.09523810 0.90476190)  
##           10) WordCount.root2< 18.9043 12   5 N (0.58333333 0.41666667) *
##           11) WordCount.root2>=18.9043 93   3 Y (0.03225806 0.96774194) *
##        3) NDSSName.my.fctrOpEd#Opinion#>=0.5 437  59 Y (0.13501144 0.86498856) *
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-50.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-51.png) 

```
##   Popular.fctr Popular.fctr.predict.Max.cor.Y.rcv.1X1.cp.0.rpart.N
## 1            N                                                3814
## 2            Y                                                 170
##   Popular.fctr.predict.Max.cor.Y.rcv.1X1.cp.0.rpart.Y
## 1                                                 127
## 2                                                 693
##          Prediction
## Reference    N    Y
##         N 3814  127
##         Y  170  693
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   9.381765e-01   7.860827e-01   9.309917e-01   9.448229e-01   8.203580e-01 
## AccuracyPValue  McnemarPValue 
##  2.798570e-127   1.480611e-02
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-52.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-53.png) 

```
##   Popular.fctr Popular.fctr.predict.Max.cor.Y.rcv.1X1.cp.0.rpart.N
## 1            N                                                1180
## 2            Y                                                  84
##   Popular.fctr.predict.Max.cor.Y.rcv.1X1.cp.0.rpart.Y
## 1                                                 318
## 2                                                 146
##          Prediction
## Reference    N    Y
##         N 1180  318
##         Y   84  146
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   7.673611e-01   2.953321e-01   7.467059e-01   7.871043e-01   8.668981e-01 
## AccuracyPValue  McnemarPValue 
##   1.000000e+00   3.224022e-31 
##                             id                            feats
## 1 Max.cor.Y.rcv.1X1.cp.0.rpart WordCount.root2,NDSSName.my.fctr
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1               0                      0.863                  0.07
##   max.AUCpROC.fit max.Sens.fit max.Spec.fit max.AUCROCR.fit
## 1       0.8821543    0.9705658    0.7937428       0.9504198
##   opt.prob.threshold.fit max.f.score.fit max.Accuracy.fit
## 1                    0.4       0.8235294        0.9381765
##   max.AccuracyLower.fit max.AccuracyUpper.fit max.Kappa.fit
## 1             0.9309917             0.9448229     0.7860827
##   max.AUCpROC.OOB max.Sens.OOB max.Spec.OOB max.AUCROCR.OOB
## 1       0.6174697    0.9218959    0.3130435       0.7773858
##   opt.prob.threshold.OOB max.f.score.OOB max.Accuracy.OOB
## 1                    0.1       0.4207493        0.7673611
##   max.AccuracyLower.OOB max.AccuracyUpper.OOB max.Kappa.OOB
## 1             0.7467059             0.7871043     0.2953321
```

```r
#stop(here"); glb_to_sav(); all.equal(glb_models_df, sav_models_df)
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
## [1] "fitting model: Max.cor.Y.rpart"
## [1] "    indep_vars: WordCount.root2,NDSSName.my.fctr"
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

![](NYTBlogs3_base_files/figure-html/fit.models_0-54.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-55.png) 

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
##              NDSSName.my.fctrOpEd#Opinion# 
##                                         55 
## NDSSName.my.fctrBusiness#Crosswords/Games# 
##                                         16 
##            NDSSName.my.fctrScience#Health# 
##                                         16 
##               NDSSName.my.fctrStyles#U.S.# 
##                                         12 
## 
## Node number 1: 4804 observations,    complexity param=0.3696408
##   predicted class=N  expected loss=0.179642  P(node) =1
##     class counts:  3941   863
##    probabilities: 0.820 0.180 
##   left son=2 (4367 obs) right son=3 (437 obs)
##   Primary splits:
##       NDSSName.my.fctrOpEd#Opinion#              < 0.5      to the left,  improve=451.59770, (0 missing)
##       NDSSName.my.fctrBusiness#Crosswords/Games# < 0.5      to the left,  improve=112.88510, (0 missing)
##       WordCount.root2                            < 25.75849 to the left,  improve=111.17610, (0 missing)
##       NDSSName.my.fctrScience#Health#            < 0.5      to the left,  improve= 99.35206, (0 missing)
##       NDSSName.my.fctrStyles#U.S.#               < 0.5      to the left,  improve= 68.73272, (0 missing)
## 
## Node number 2: 4367 observations,    complexity param=0.09849363
##   predicted class=N  expected loss=0.1110602  P(node) =0.9090341
##     class counts:  3882   485
##    probabilities: 0.889 0.111 
##   left son=4 (4262 obs) right son=5 (105 obs)
##   Primary splits:
##       NDSSName.my.fctrBusiness#Crosswords/Games# < 0.5      to the left,  improve=135.55130, (0 missing)
##       NDSSName.my.fctrScience#Health#            < 0.5      to the left,  improve=125.07920, (0 missing)
##       WordCount.root2                            < 25.75849 to the left,  improve= 94.70710, (0 missing)
##       NDSSName.my.fctrStyles#U.S.#               < 0.5      to the left,  improve= 88.56821, (0 missing)
##       NDSSName.my.fctr#Opinion#ThePublicEditor   < 0.5      to the left,  improve= 18.74400, (0 missing)
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
##       NDSSName.my.fctrScience#Health#          < 0.5      to the left,  improve=132.96710, (0 missing)
##       NDSSName.my.fctrStyles#U.S.#             < 0.5      to the left,  improve= 94.69099, (0 missing)
##       WordCount.root2                          < 26.49528 to the left,  improve= 84.07487, (0 missing)
##       NDSSName.my.fctr#Opinion#ThePublicEditor < 0.5      to the left,  improve= 19.71762, (0 missing)
##       NDSSName.my.fctrTStyle##                 < 0.5      to the right, improve= 10.17000, (0 missing)
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
##       NDSSName.my.fctrStyles#U.S.#             < 0.5      to the left,  improve=102.410700, (0 missing)
##       WordCount.root2                          < 25.01    to the left,  improve= 47.352210, (0 missing)
##       NDSSName.my.fctr#Opinion#ThePublicEditor < 0.5      to the left,  improve= 20.930810, (0 missing)
##       NDSSName.my.fctrTStyle##                 < 0.5      to the right, improve=  5.249425, (0 missing)
##       NDSSName.my.fctrBusiness#Technology#     < 0.5      to the left,  improve=  2.395935, (0 missing)
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
##    2) NDSSName.my.fctrOpEd#Opinion#< 0.5 4367 485 N (0.88893978 0.11106022)  
##      4) NDSSName.my.fctrBusiness#Crosswords/Games#< 0.5 4262 390 N (0.90849366 0.09150634)  
##        8) NDSSName.my.fctrScience#Health#< 0.5 4114 279 N (0.93218279 0.06781721)  
##         16) NDSSName.my.fctrStyles#U.S.#< 0.5 3987 191 N (0.95209431 0.04790569) *
##         17) NDSSName.my.fctrStyles#U.S.#>=0.5 127  39 Y (0.30708661 0.69291339) *
##        9) NDSSName.my.fctrScience#Health#>=0.5 148  37 Y (0.25000000 0.75000000) *
##      5) NDSSName.my.fctrBusiness#Crosswords/Games#>=0.5 105  10 Y (0.09523810 0.90476190) *
##    3) NDSSName.my.fctrOpEd#Opinion#>=0.5 437  59 Y (0.13501144 0.86498856) *
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-56.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-57.png) 

```
##   Popular.fctr Popular.fctr.predict.Max.cor.Y.rpart.N
## 1            N                                   3796
## 2            Y                                    191
##   Popular.fctr.predict.Max.cor.Y.rpart.Y
## 1                                    145
## 2                                    672
##          Prediction
## Reference    N    Y
##         N 3796  145
##         Y  191  672
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   9.300583e-01   7.576571e-01   9.224771e-01   9.371115e-01   8.203580e-01 
## AccuracyPValue  McnemarPValue 
##  4.458834e-108   1.409037e-02
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-58.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-59.png) 

```
##   Popular.fctr Popular.fctr.predict.Max.cor.Y.rpart.N
## 1            N                                   1355
## 2            Y                                    168
##   Popular.fctr.predict.Max.cor.Y.rpart.Y
## 1                                    143
## 2                                     62
##          Prediction
## Reference    N    Y
##         N 1355  143
##         Y  168   62
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##      0.8200231      0.1825002      0.8010821      0.8378705      0.8668981 
## AccuracyPValue  McnemarPValue 
##      1.0000000      0.1735405 
##                id                            feats max.nTuningRuns
## 1 Max.cor.Y.rpart WordCount.root2,NDSSName.my.fctr               5
##   min.elapsedtime.everything min.elapsedtime.final max.AUCpROC.fit
## 1                      2.882                 0.071       0.8709432
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
if (!is.null(glb_date_vars) && 
    (sum(grepl(paste(glb_date_vars, "\\.day\\.minutes\\.poly\\.", sep=""),
               names(glbObsAll))) > 0)) {
# ret_lst <- myfit_mdl(mdl_id="Max.cor.Y.TmSrs.poly1", 
#                         model_method=ifelse(glb_is_regression, "lm", 
#                                         ifelse(glb_is_binomial, "glm", "rpart")),
#                      model_type=glb_model_type,
#                         indep_vars_vctr=c(max_cor_y_x_vars, paste0(glb_date_vars, ".day.minutes")),
#                         rsp_var=glb_rsp_var, rsp_var_out=glb_rsp_var_out,
#                         fit_df=glbObsFit, OOB_df=glbObsOOB,
#                         n_cv_folds=glb_rcv_n_folds, tune_models_df=NULL)
# 
ret_lst <- myfit_mdl(mdl_id="Max.cor.Y.TmSrs.poly", 
                        model_method=ifelse(glb_is_regression, "lm", 
                                        ifelse(glb_is_binomial, "glm", "rpart")),
                     model_type=glb_model_type,
                        indep_vars_vctr=c(max_cor_y_x_vars, 
            grep(paste(glb_date_vars, "\\.day\\.minutes\\.poly\\.", sep=""),
                        names(glbObsAll), value=TRUE)),
                        rsp_var=glb_rsp_var, rsp_var_out=glb_rsp_var_out,
                        fit_df=glbObsFit, OOB_df=glbObsOOB,
                        n_cv_folds=glb_rcv_n_folds, tune_models_df=NULL)
}

# Interactions.High.cor.Y
if (length(int_feats <- setdiff(setdiff(unique(glb_feats_df$cor.high.X), NA), 
                                subset(glb_feats_df, nzv)$id)) > 0) {
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
## [1] "fitting model: Interact.High.cor.Y.glmnet"
## [1] "    indep_vars: WordCount.root2,NDSSName.my.fctr,WordCount.root2:WordCount.root2"
## Aggregating results
## Selecting tuning parameters
## Fitting alpha = 0.325, lambda = 0.0201 on full training set
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-60.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-61.png) 

```
##             Length Class      Mode     
## a0            99   -none-     numeric  
## beta        2079   dgCMatrix  S4       
## df            99   -none-     numeric  
## dim            2   -none-     numeric  
## lambda        99   -none-     numeric  
## dev.ratio     99   -none-     numeric  
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
##                                        (Intercept) 
##                                        -3.89350373 
##                       NDSSName.my.fctr#Multimedia# 
##                                        -0.01916344 
##             NDSSName.my.fctr#Opinion#RoomForDebate 
##                                        -1.18453357 
##           NDSSName.my.fctr#Opinion#ThePublicEditor 
##                                         3.21701058 
##                    NDSSName.my.fctr#U.S.#Education 
##                                        -0.47679040 
## NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness 
##                                        -0.09891374 
##         NDSSName.my.fctrBusiness#Crosswords/Games# 
##                                         3.87281404 
##               NDSSName.my.fctrBusiness#Technology# 
##                                         0.40965256 
##                     NDSSName.my.fctrForeign#World# 
##                                        -0.05114617 
##          NDSSName.my.fctrForeign#World#AsiaPacific 
##                                        -0.47464340 
##                      NDSSName.my.fctrOpEd#Opinion# 
##                                         3.95214357 
##                    NDSSName.my.fctrScience#Health# 
##                                         3.14232408 
##                    NDSSName.my.fctrStyles##Fashion 
##                                        -0.31867093 
##                       NDSSName.my.fctrStyles#U.S.# 
##                                         2.92610567 
##                           NDSSName.my.fctrTStyle## 
##                                        -0.60538025 
##                                    WordCount.root2 
##                                         0.05783392 
## [1] "max lambda < lambdaOpt:"
##                                        (Intercept) 
##                                        -3.95644632 
##                       NDSSName.my.fctr#Multimedia# 
##                                        -0.07182859 
##             NDSSName.my.fctr#Opinion#RoomForDebate 
##                                        -1.30034382 
##           NDSSName.my.fctr#Opinion#ThePublicEditor 
##                                         3.29694236 
##                    NDSSName.my.fctr#U.S.#Education 
##                                        -0.53415905 
## NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness 
##                                        -0.14259759 
##         NDSSName.my.fctrBusiness#Crosswords/Games# 
##                                         3.94231812 
##               NDSSName.my.fctrBusiness#Technology# 
##                                         0.45657914 
##                     NDSSName.my.fctrForeign#World# 
##                                        -0.10021084 
##          NDSSName.my.fctrForeign#World#AsiaPacific 
##                                        -0.53077048 
##                      NDSSName.my.fctrOpEd#Opinion# 
##                                         4.01324666 
##                    NDSSName.my.fctrScience#Health# 
##                                         3.18936803 
##                    NDSSName.my.fctrStyles##Fashion 
##                                        -0.38069674 
##                       NDSSName.my.fctrStyles#U.S.# 
##                                         2.97176051 
##                           NDSSName.my.fctrTStyle## 
##                                        -0.64837249 
##                                    WordCount.root2 
##                                         0.05978318
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-62.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-63.png) 

```
##   Popular.fctr Popular.fctr.predict.Interact.High.cor.Y.glmnet.N
## 1            N                                              3796
## 2            Y                                               177
##   Popular.fctr.predict.Interact.High.cor.Y.glmnet.Y
## 1                                               145
## 2                                               686
##          Prediction
## Reference    N    Y
##         N 3796  145
##         Y  177  686
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   9.329725e-01   7.692476e-01   9.255302e-01   9.398832e-01   8.203580e-01 
## AccuracyPValue  McnemarPValue 
##  1.026390e-114   8.406670e-02
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-64.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-65.png) 

```
##   Popular.fctr Popular.fctr.predict.Interact.High.cor.Y.glmnet.N
## 1            N                                              1146
## 2            Y                                                67
##   Popular.fctr.predict.Interact.High.cor.Y.glmnet.Y
## 1                                               352
## 2                                               163
##          Prediction
## Reference    N    Y
##         N 1146  352
##         Y   67  163
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   7.575231e-01   3.107477e-01   7.365992e-01   7.775689e-01   8.668981e-01 
## AccuracyPValue  McnemarPValue 
##   1.000000e+00   9.066396e-44 
##                           id
## 1 Interact.High.cor.Y.glmnet
##                                                              feats
## 1 WordCount.root2,NDSSName.my.fctr,WordCount.root2:WordCount.root2
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1              25                      4.395                  0.27
##   max.AUCpROC.fit max.Sens.fit max.Spec.fit max.AUCROCR.fit
## 1       0.8767919     0.964476    0.7891078       0.9582555
##   opt.prob.threshold.fit max.f.score.fit max.Accuracy.fit
## 1                    0.4       0.8099174        0.9333193
##   max.AccuracyLower.fit max.AccuracyUpper.fit max.Kappa.fit
## 1             0.9255302             0.9398832     0.7690803
##   max.AUCpROC.OOB max.Sens.OOB max.Spec.OOB max.AUCROCR.OOB
## 1       0.5962443    0.9098798    0.2826087       0.8067975
##   opt.prob.threshold.OOB max.f.score.OOB max.Accuracy.OOB
## 1                    0.1       0.4375839        0.7575231
##   max.AccuracyLower.OOB max.AccuracyUpper.OOB max.Kappa.OOB
## 1             0.7365992             0.7775689     0.3107477
##   max.AccuracySD.fit max.KappaSD.fit
## 1        0.005178375      0.01754365
```

```r
# Low.cor.X
# if (glb_is_classification && glb_is_binomial)
#     indep_vars_vctr <- subset(glb_feats_df, is.na(cor.high.X) & 
#                                             is.ConditionalX.y & 
#                                             (exclude.as.feat != 1))[, "id"] else
indep_vars <- subset(glb_feats_df, is.na(cor.high.X) & !nzv & 
                              (exclude.as.feat != 1))[, "id"]  
indep_vars <- myadjust_interaction_feats(indep_vars)
ret_lst <- myfit_mdl(mdl_specs_lst=myinit_mdl_specs_lst(mdl_specs_lst=list(
        id.prefix="Low.cor.X", 
        type=glb_model_type, trainControl.method="repeatedcv",
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
## [1] "fitting model: Low.cor.X.glmnet"
## [1] "    indep_vars: WordCount.root2,NDSSName.my.fctr,.rnorm,WordCount.nexp"
## Aggregating results
## Selecting tuning parameters
## Fitting alpha = 0.325, lambda = 0.0201 on full training set
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-66.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-67.png) 

```
##             Length Class      Mode     
## a0           100   -none-     numeric  
## beta        2300   dgCMatrix  S4       
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
## xNames        23   -none-     character
## problemType    1   -none-     character
## tuneValue      2   data.frame list     
## obsLevels      2   -none-     character
## [1] "min lambda > lambdaOpt:"
##                                        (Intercept) 
##                                        -3.89350373 
##                       NDSSName.my.fctr#Multimedia# 
##                                        -0.01916344 
##             NDSSName.my.fctr#Opinion#RoomForDebate 
##                                        -1.18453357 
##           NDSSName.my.fctr#Opinion#ThePublicEditor 
##                                         3.21701058 
##                    NDSSName.my.fctr#U.S.#Education 
##                                        -0.47679040 
## NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness 
##                                        -0.09891374 
##         NDSSName.my.fctrBusiness#Crosswords/Games# 
##                                         3.87281404 
##               NDSSName.my.fctrBusiness#Technology# 
##                                         0.40965256 
##                     NDSSName.my.fctrForeign#World# 
##                                        -0.05114617 
##          NDSSName.my.fctrForeign#World#AsiaPacific 
##                                        -0.47464340 
##                      NDSSName.my.fctrOpEd#Opinion# 
##                                         3.95214357 
##                    NDSSName.my.fctrScience#Health# 
##                                         3.14232408 
##                    NDSSName.my.fctrStyles##Fashion 
##                                        -0.31867093 
##                       NDSSName.my.fctrStyles#U.S.# 
##                                         2.92610567 
##                           NDSSName.my.fctrTStyle## 
##                                        -0.60538025 
##                                    WordCount.root2 
##                                         0.05783392 
## [1] "max lambda < lambdaOpt:"
##                                        (Intercept) 
##                                        -3.95644632 
##                       NDSSName.my.fctr#Multimedia# 
##                                        -0.07182859 
##             NDSSName.my.fctr#Opinion#RoomForDebate 
##                                        -1.30034382 
##           NDSSName.my.fctr#Opinion#ThePublicEditor 
##                                         3.29694236 
##                    NDSSName.my.fctr#U.S.#Education 
##                                        -0.53415905 
## NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness 
##                                        -0.14259759 
##         NDSSName.my.fctrBusiness#Crosswords/Games# 
##                                         3.94231812 
##               NDSSName.my.fctrBusiness#Technology# 
##                                         0.45657914 
##                     NDSSName.my.fctrForeign#World# 
##                                        -0.10021084 
##          NDSSName.my.fctrForeign#World#AsiaPacific 
##                                        -0.53077048 
##                      NDSSName.my.fctrOpEd#Opinion# 
##                                         4.01324666 
##                    NDSSName.my.fctrScience#Health# 
##                                         3.18936803 
##                    NDSSName.my.fctrStyles##Fashion 
##                                        -0.38069674 
##                       NDSSName.my.fctrStyles#U.S.# 
##                                         2.97176051 
##                           NDSSName.my.fctrTStyle## 
##                                        -0.64837249 
##                                    WordCount.root2 
##                                         0.05978318
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-68.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-69.png) 

```
##   Popular.fctr Popular.fctr.predict.Low.cor.X.glmnet.N
## 1            N                                    3796
## 2            Y                                     177
##   Popular.fctr.predict.Low.cor.X.glmnet.Y
## 1                                     145
## 2                                     686
##          Prediction
## Reference    N    Y
##         N 3796  145
##         Y  177  686
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   9.329725e-01   7.692476e-01   9.255302e-01   9.398832e-01   8.203580e-01 
## AccuracyPValue  McnemarPValue 
##  1.026390e-114   8.406670e-02
```

![](NYTBlogs3_base_files/figure-html/fit.models_0-70.png) ![](NYTBlogs3_base_files/figure-html/fit.models_0-71.png) 

```
##   Popular.fctr Popular.fctr.predict.Low.cor.X.glmnet.N
## 1            N                                    1146
## 2            Y                                      67
##   Popular.fctr.predict.Low.cor.X.glmnet.Y
## 1                                     352
## 2                                     163
##          Prediction
## Reference    N    Y
##         N 1146  352
##         Y   67  163
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   7.575231e-01   3.107477e-01   7.365992e-01   7.775689e-01   8.668981e-01 
## AccuracyPValue  McnemarPValue 
##   1.000000e+00   9.066396e-44 
##                 id                                                  feats
## 1 Low.cor.X.glmnet WordCount.root2,NDSSName.my.fctr,.rnorm,WordCount.nexp
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1              25                      4.549                 0.293
##   max.AUCpROC.fit max.Sens.fit max.Spec.fit max.AUCROCR.fit
## 1       0.8767919     0.964476    0.7891078       0.9582555
##   opt.prob.threshold.fit max.f.score.fit max.Accuracy.fit
## 1                    0.4       0.8099174        0.9333193
##   max.AccuracyLower.fit max.AccuracyUpper.fit max.Kappa.fit
## 1             0.9255302             0.9398832     0.7690803
##   max.AUCpROC.OOB max.Sens.OOB max.Spec.OOB max.AUCROCR.OOB
## 1       0.5962443    0.9098798    0.2826087       0.8067975
##   opt.prob.threshold.OOB max.f.score.OOB max.Accuracy.OOB
## 1                    0.1       0.4375839        0.7575231
##   max.AccuracyLower.OOB max.AccuracyUpper.OOB max.Kappa.OOB
## 1             0.7365992             0.7775689     0.3107477
##   max.AccuracySD.fit max.KappaSD.fit
## 1        0.005178375      0.01754365
```

```r
rm(ret_lst)

glb_chunks_df <- myadd_chunk(glb_chunks_df, "fit.models", major.inc=FALSE)
```

```
##         label step_major step_minor label_minor     bgn    end elapsed
## 10 fit.models          6          0           0  35.238 139.59 104.352
## 11 fit.models          6          1           1 139.591     NA      NA
```


```r
fit.models_1_chunk_df <- myadd_chunk(NULL, "fit.models_1_bgn", label.minor="setup")
```

```
##              label step_major step_minor label_minor     bgn end elapsed
## 1 fit.models_1_bgn          1          0       setup 151.043  NA      NA
```

```r
#stop(here"); glb_to_sav(); all.equal(glb_models_df, sav_models_df)
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
        #   select most importance feature
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
            importance_cutoff <- max(10, imp_rnorm, na.rm=TRUE)
            interact_vars <- 
                tail(row.names(subset(bst_featsimp_df, 
                                      importance > importance_cutoff)), -1)
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
    
    if (grepl("^%<d-%", indep_vars)) {    
#stop(here")        
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
        ret_lst <- 
            myfit_mdl(mdl_specs_lst = myinit_mdl_specs_lst(mdl_specs_lst = list(
            id.prefix = mdl_id_pfx, 
            type = glb_model_type, tune.df = glb_tune_models_df,
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
## 1   fit.models_1_bgn          1          0       setup 151.043 151.052
## 2 fit.models_1_All.X          1          1       setup 151.052      NA
##   elapsed
## 1   0.009
## 2      NA
```

```
## Warning in if (grepl("^%<d-%", indep_vars)) {: the condition has length > 1
## and only the first element will be used
```

```
##                label step_major step_minor label_minor     bgn     end
## 2 fit.models_1_All.X          1          1       setup 151.052 151.059
## 3 fit.models_1_All.X          1          2      glmnet 151.060      NA
##   elapsed
## 2   0.007
## 3      NA
## [1] "fitting model: All.X.glmnet"
## [1] "    indep_vars: WordCount.root2,WordCount.log1p,NDSSName.my.fctr,.rnorm,WordCount.nexp"
## Aggregating results
## Selecting tuning parameters
## Fitting alpha = 1, lambda = 0.00434 on full training set
```

```
## Warning in myfit_mdl(mdl_specs_lst = myinit_mdl_specs_lst(mdl_specs_lst
## = list(id.prefix = mdl_id_pfx, : model's bestTune found at an extreme of
## tuneGrid for parameter: alpha
```

![](NYTBlogs3_base_files/figure-html/fit.models_1-1.png) ![](NYTBlogs3_base_files/figure-html/fit.models_1-2.png) 

```
##             Length Class      Mode     
## a0            91   -none-     numeric  
## beta        2184   dgCMatrix  S4       
## df            91   -none-     numeric  
## dim            2   -none-     numeric  
## lambda        91   -none-     numeric  
## dev.ratio     91   -none-     numeric  
## nulldev        1   -none-     numeric  
## npasses        1   -none-     numeric  
## jerr           1   -none-     numeric  
## offset         1   -none-     logical  
## classnames     2   -none-     character
## call           5   -none-     call     
## nobs           1   -none-     numeric  
## lambdaOpt      1   -none-     numeric  
## xNames        24   -none-     character
## problemType    1   -none-     character
## tuneValue      2   data.frame list     
## obsLevels      2   -none-     character
## [1] "min lambda > lambdaOpt:"
##                                        (Intercept) 
##                                        -5.52332827 
##                       NDSSName.my.fctr#Multimedia# 
##                                        -0.16474095 
##             NDSSName.my.fctr#Opinion#RoomForDebate 
##                                        -2.18132553 
##           NDSSName.my.fctr#Opinion#ThePublicEditor 
##                                         3.90787932 
##                    NDSSName.my.fctr#U.S.#Education 
##                                        -0.90764328 
## NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness 
##                                        -0.24502067 
##         NDSSName.my.fctrBusiness#Crosswords/Games# 
##                                         4.67941504 
##               NDSSName.my.fctrBusiness#Technology# 
##                                         0.83419892 
##                     NDSSName.my.fctrForeign#World# 
##                                        -0.08805297 
##          NDSSName.my.fctrForeign#World#AsiaPacific 
##                                        -0.83659624 
##                      NDSSName.my.fctrOpEd#Opinion# 
##                                         4.79265681 
##                    NDSSName.my.fctrScience#Health# 
##                                         3.67959505 
##                    NDSSName.my.fctrStyles##Fashion 
##                                        -0.64461547 
##                       NDSSName.my.fctrStyles#U.S.# 
##                                         3.44343706 
##                           NDSSName.my.fctrTStyle## 
##                                        -0.97236910 
##                            NDSSName.my.fctrmyOther 
##                                        -0.23762498 
##                                    WordCount.log1p 
##                                         0.18619880 
##                                    WordCount.root2 
##                                         0.06691104 
## [1] "max lambda < lambdaOpt:"
##                                        (Intercept) 
##                                        -5.59881022 
##                       NDSSName.my.fctr#Multimedia# 
##                                        -0.24609678 
##             NDSSName.my.fctr#Opinion#RoomForDebate 
##                                        -2.30818954 
##           NDSSName.my.fctr#Opinion#ThePublicEditor 
##                                         3.94273791 
##                    NDSSName.my.fctr#U.S.#Education 
##                                        -1.01466873 
## NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness 
##                                        -0.30479819 
##         NDSSName.my.fctrBusiness#Crosswords/Games# 
##                                         4.69839363 
##               NDSSName.my.fctrBusiness#Technology# 
##                                         0.84455935 
##                     NDSSName.my.fctrForeign#World# 
##                                        -0.18204983 
##          NDSSName.my.fctrForeign#World#AsiaPacific 
##                                        -0.91223472 
##                      NDSSName.my.fctrOpEd#Opinion# 
##                                         4.80649138 
##                    NDSSName.my.fctrScience#Health# 
##                                         3.68561250 
##                    NDSSName.my.fctrStyles##Fashion 
##                                        -0.74721366 
##                       NDSSName.my.fctrStyles#U.S.# 
##                                         3.44963924 
##                           NDSSName.my.fctrTStyle## 
##                                        -1.04099014 
##                     NDSSName.my.fctrTravel#Travel# 
##                                        -0.03203556 
##                            NDSSName.my.fctrmyOther 
##                                        -0.35136209 
##                                    WordCount.log1p 
##                                         0.19889748 
##                                    WordCount.root2 
##                                         0.06700663
```

![](NYTBlogs3_base_files/figure-html/fit.models_1-3.png) ![](NYTBlogs3_base_files/figure-html/fit.models_1-4.png) 

```
##   Popular.fctr Popular.fctr.predict.All.X.glmnet.N
## 1            N                                3793
## 2            Y                                 176
##   Popular.fctr.predict.All.X.glmnet.Y
## 1                                 148
## 2                                 687
##          Prediction
## Reference    N    Y
##         N 3793  148
##         Y  176  687
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   9.325562e-01   7.682400e-01   9.250938e-01   9.394875e-01   8.203580e-01 
## AccuracyPValue  McnemarPValue 
##  9.476102e-114   1.336144e-01
```

![](NYTBlogs3_base_files/figure-html/fit.models_1-5.png) ![](NYTBlogs3_base_files/figure-html/fit.models_1-6.png) 

```
##   Popular.fctr Popular.fctr.predict.All.X.glmnet.N
## 1            N                                1176
## 2            Y                                  77
##   Popular.fctr.predict.All.X.glmnet.Y
## 1                                 322
## 2                                 153
##          Prediction
## Reference    N    Y
##         N 1176  322
##         Y   77  153
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   7.690972e-01   3.103487e-01   7.484910e-01   7.887855e-01   8.668981e-01 
## AccuracyPValue  McnemarPValue 
##   1.000000e+00   2.576173e-34 
##             id
## 1 All.X.glmnet
##                                                                    feats
## 1 WordCount.root2,WordCount.log1p,NDSSName.my.fctr,.rnorm,WordCount.nexp
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1              25                      4.817                 0.294
##   max.AUCpROC.fit max.Sens.fit max.Spec.fit max.AUCROCR.fit
## 1       0.8771175    0.9639685    0.7902665       0.9589072
##   opt.prob.threshold.fit max.f.score.fit max.Accuracy.fit
## 1                    0.4       0.8091873        0.9326952
##   max.AccuracyLower.fit max.AccuracyUpper.fit max.Kappa.fit
## 1             0.9250938             0.9394875     0.7668828
##   max.AUCpROC.OOB max.Sens.OOB max.Spec.OOB max.AUCROCR.OOB
## 1        0.596578    0.9105474    0.2826087       0.8061009
##   opt.prob.threshold.OOB max.f.score.OOB max.Accuracy.OOB
## 1                    0.1       0.4340426        0.7690972
##   max.AccuracyLower.OOB max.AccuracyUpper.OOB max.Kappa.OOB
## 1              0.748491             0.7887855     0.3103487
##   max.AccuracySD.fit max.KappaSD.fit
## 1        0.004892081      0.01659075
```

```r
# Check if other preProcess methods improve model performance
fit.models_1_chunk_df <- 
    myadd_chunk(fit.models_1_chunk_df, "fit.models_1_preProc", major.inc = FALSE,
                label.minor = "preProc")
```

```
##                  label step_major step_minor label_minor     bgn     end
## 3   fit.models_1_All.X          1          2      glmnet 151.060 161.253
## 4 fit.models_1_preProc          1          3     preProc 161.253      NA
##   elapsed
## 3  10.193
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
            type=glb_model_type, tune.df=glb_tune_models_df,
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
#                                 rsp_var=glb_rsp_var, rsp_var_out=glb_rsp_var_out,
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
#                                 rsp_var=glb_rsp_var, rsp_var_out=glb_rsp_var_out,
#                                 fit_df=glbObsFit, OOB_df=glbObsOOB,
#                     n_cv_folds=glb_rcv_n_folds, tune_models_df=glb_tune_models_df)
#     csm_mdl_id <- paste0(mdl_id, ".", method)
#     csm_featsimp_df <- myget_feats_importance(glb_models_lst[[paste0(mdl_id, ".",
#                                                                      method)]]);               print(head(csm_featsimp_df))
# }
###

# Ntv.1.lm <- lm(reformulate(indep_vars_vctr, glb_rsp_var), glbObsTrn); print(summary(Ntv.1.lm))

#glb_models_df[, "max.Accuracy.OOB", FALSE]
#varImp(glb_models_lst[["Low.cor.X.glm"]])
#orderBy(~ -Overall, varImp(glb_models_lst[["All.X.2.glm"]])$importance)
#orderBy(~ -Overall, varImp(glb_models_lst[["All.X.3.glm"]])$importance)
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
#                             rsp_var=glb_rsp_var, rsp_var_out=glb_rsp_var_out,
#                             fit_df=glbObsFit, OOB_df=glbObsOOB,
#                             n_cv_folds=glb_rcv_n_folds, tune_models_df=glb_tune_models_df,
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
##                                                        id
## MFO.myMFO_classfr                       MFO.myMFO_classfr
## Random.myrandom_classfr           Random.myrandom_classfr
## Max.cor.Y.rcv.1X1.glmnet         Max.cor.Y.rcv.1X1.glmnet
## Max.cor.Y.rcv.3X1.glmnet         Max.cor.Y.rcv.3X1.glmnet
## Max.cor.Y.rcv.3X3.glmnet         Max.cor.Y.rcv.3X3.glmnet
## Max.cor.Y.rcv.3X5.glmnet         Max.cor.Y.rcv.3X5.glmnet
## Max.cor.Y.rcv.5X1.glmnet         Max.cor.Y.rcv.5X1.glmnet
## Max.cor.Y.rcv.5X3.glmnet         Max.cor.Y.rcv.5X3.glmnet
## Max.cor.Y.rcv.5X5.glmnet         Max.cor.Y.rcv.5X5.glmnet
## Max.cor.Y.rcv.1X1.cp.0.rpart Max.cor.Y.rcv.1X1.cp.0.rpart
## Max.cor.Y.rpart                           Max.cor.Y.rpart
## Interact.High.cor.Y.glmnet     Interact.High.cor.Y.glmnet
## Low.cor.X.glmnet                         Low.cor.X.glmnet
## All.X.glmnet                                 All.X.glmnet
##                                                                                               feats
## MFO.myMFO_classfr                                                                            .rnorm
## Random.myrandom_classfr                                                                      .rnorm
## Max.cor.Y.rcv.1X1.glmnet                                           WordCount.root2,NDSSName.my.fctr
## Max.cor.Y.rcv.3X1.glmnet                                           WordCount.root2,NDSSName.my.fctr
## Max.cor.Y.rcv.3X3.glmnet                                           WordCount.root2,NDSSName.my.fctr
## Max.cor.Y.rcv.3X5.glmnet                                           WordCount.root2,NDSSName.my.fctr
## Max.cor.Y.rcv.5X1.glmnet                                           WordCount.root2,NDSSName.my.fctr
## Max.cor.Y.rcv.5X3.glmnet                                           WordCount.root2,NDSSName.my.fctr
## Max.cor.Y.rcv.5X5.glmnet                                           WordCount.root2,NDSSName.my.fctr
## Max.cor.Y.rcv.1X1.cp.0.rpart                                       WordCount.root2,NDSSName.my.fctr
## Max.cor.Y.rpart                                                    WordCount.root2,NDSSName.my.fctr
## Interact.High.cor.Y.glmnet         WordCount.root2,NDSSName.my.fctr,WordCount.root2:WordCount.root2
## Low.cor.X.glmnet                             WordCount.root2,NDSSName.my.fctr,.rnorm,WordCount.nexp
## All.X.glmnet                 WordCount.root2,WordCount.log1p,NDSSName.my.fctr,.rnorm,WordCount.nexp
##                              max.nTuningRuns min.elapsedtime.everything
## MFO.myMFO_classfr                          0                      0.269
## Random.myrandom_classfr                    0                      0.269
## Max.cor.Y.rcv.1X1.glmnet                   0                      0.972
## Max.cor.Y.rcv.3X1.glmnet                  25                      2.826
## Max.cor.Y.rcv.3X3.glmnet                  25                      4.817
## Max.cor.Y.rcv.3X5.glmnet                  25                      6.419
## Max.cor.Y.rcv.5X1.glmnet                  25                      3.223
## Max.cor.Y.rcv.5X3.glmnet                  25                      6.142
## Max.cor.Y.rcv.5X5.glmnet                  25                      8.908
## Max.cor.Y.rcv.1X1.cp.0.rpart               0                      0.863
## Max.cor.Y.rpart                            5                      2.882
## Interact.High.cor.Y.glmnet                25                      4.395
## Low.cor.X.glmnet                          25                      4.549
## All.X.glmnet                              25                      4.817
##                              min.elapsedtime.final max.AUCpROC.fit
## MFO.myMFO_classfr                            0.003       0.5000000
## Random.myrandom_classfr                      0.001       0.4990604
## Max.cor.Y.rcv.1X1.glmnet                     0.273       0.8790544
## Max.cor.Y.rcv.3X1.glmnet                     0.270       0.8767919
## Max.cor.Y.rcv.3X3.glmnet                     0.270       0.8767919
## Max.cor.Y.rcv.3X5.glmnet                     0.269       0.8767919
## Max.cor.Y.rcv.5X1.glmnet                     0.267       0.8784031
## Max.cor.Y.rcv.5X3.glmnet                     0.267       0.8784031
## Max.cor.Y.rcv.5X5.glmnet                     0.268       0.8784031
## Max.cor.Y.rcv.1X1.cp.0.rpart                 0.070       0.8821543
## Max.cor.Y.rpart                              0.071       0.8709432
## Interact.High.cor.Y.glmnet                   0.270       0.8767919
## Low.cor.X.glmnet                             0.293       0.8767919
## All.X.glmnet                                 0.294       0.8771175
##                              max.Sens.fit max.Spec.fit max.AUCROCR.fit
## MFO.myMFO_classfr               1.0000000    0.0000000       0.5000000
## Random.myrandom_classfr         0.8312611    0.1668598       0.4972757
## Max.cor.Y.rcv.1X1.glmnet        0.9632073    0.7949015       0.9608594
## Max.cor.Y.rcv.3X1.glmnet        0.9644760    0.7891078       0.9582555
## Max.cor.Y.rcv.3X3.glmnet        0.9644760    0.7891078       0.9582555
## Max.cor.Y.rcv.3X5.glmnet        0.9644760    0.7891078       0.9582555
## Max.cor.Y.rcv.5X1.glmnet        0.9642223    0.7925840       0.9607052
## Max.cor.Y.rcv.5X3.glmnet        0.9642223    0.7925840       0.9607052
## Max.cor.Y.rcv.5X5.glmnet        0.9642223    0.7925840       0.9607052
## Max.cor.Y.rcv.1X1.cp.0.rpart    0.9705658    0.7937428       0.9504198
## Max.cor.Y.rpart                 0.9632073    0.7786790       0.8746354
## Interact.High.cor.Y.glmnet      0.9644760    0.7891078       0.9582555
## Low.cor.X.glmnet                0.9644760    0.7891078       0.9582555
## All.X.glmnet                    0.9639685    0.7902665       0.9589072
##                              opt.prob.threshold.fit max.f.score.fit
## MFO.myMFO_classfr                               0.1       0.3045703
## Random.myrandom_classfr                         0.1       0.3045703
## Max.cor.Y.rcv.1X1.glmnet                        0.5       0.8099174
## Max.cor.Y.rcv.3X1.glmnet                        0.4       0.8099174
## Max.cor.Y.rcv.3X3.glmnet                        0.4       0.8099174
## Max.cor.Y.rcv.3X5.glmnet                        0.4       0.8099174
## Max.cor.Y.rcv.5X1.glmnet                        0.5       0.8104265
## Max.cor.Y.rcv.5X3.glmnet                        0.5       0.8104265
## Max.cor.Y.rcv.5X5.glmnet                        0.5       0.8104265
## Max.cor.Y.rcv.1X1.cp.0.rpart                    0.4       0.8235294
## Max.cor.Y.rpart                                 0.6       0.8000000
## Interact.High.cor.Y.glmnet                      0.4       0.8099174
## Low.cor.X.glmnet                                0.4       0.8099174
## All.X.glmnet                                    0.4       0.8091873
##                              max.Accuracy.fit max.AccuracyLower.fit
## MFO.myMFO_classfr                   0.1796420             0.1688795
## Random.myrandom_classfr             0.1796420             0.1688795
## Max.cor.Y.rcv.1X1.glmnet            0.9329725             0.9255302
## Max.cor.Y.rcv.3X1.glmnet            0.9335973             0.9255302
## Max.cor.Y.rcv.3X3.glmnet            0.9333193             0.9255302
## Max.cor.Y.rcv.3X5.glmnet            0.9332218             0.9255302
## Max.cor.Y.rcv.5X1.glmnet            0.9331818             0.9259666
## Max.cor.Y.rcv.5X3.glmnet            0.9333905             0.9259666
## Max.cor.Y.rcv.5X5.glmnet            0.9331816             0.9259666
## Max.cor.Y.rcv.1X1.cp.0.rpart        0.9381765             0.9309917
## Max.cor.Y.rpart                     0.9296422             0.9224771
## Interact.High.cor.Y.glmnet          0.9333193             0.9255302
## Low.cor.X.glmnet                    0.9333193             0.9255302
## All.X.glmnet                        0.9326952             0.9250938
##                              max.AccuracyUpper.fit max.Kappa.fit
## MFO.myMFO_classfr                        0.1907952     0.0000000
## Random.myrandom_classfr                  0.1907952     0.0000000
## Max.cor.Y.rcv.1X1.glmnet                 0.9398832     0.7692476
## Max.cor.Y.rcv.3X1.glmnet                 0.9398832     0.7691678
## Max.cor.Y.rcv.3X3.glmnet                 0.9398832     0.7690803
## Max.cor.Y.rcv.3X5.glmnet                 0.9398832     0.7686375
## Max.cor.Y.rcv.5X1.glmnet                 0.9402789     0.7689055
## Max.cor.Y.rcv.5X3.glmnet                 0.9402789     0.7698577
## Max.cor.Y.rcv.5X5.glmnet                 0.9402789     0.7691429
## Max.cor.Y.rcv.1X1.cp.0.rpart             0.9448229     0.7860827
## Max.cor.Y.rpart                          0.9371115     0.7515134
## Interact.High.cor.Y.glmnet               0.9398832     0.7690803
## Low.cor.X.glmnet                         0.9398832     0.7690803
## All.X.glmnet                             0.9394875     0.7668828
##                              max.AUCpROC.OOB max.Sens.OOB max.Spec.OOB
## MFO.myMFO_classfr                  0.5000000    1.0000000    0.0000000
## Random.myrandom_classfr            0.5125675    0.8077437    0.2173913
## Max.cor.Y.rcv.1X1.glmnet           0.5962443    0.9098798    0.2826087
## Max.cor.Y.rcv.3X1.glmnet           0.5962443    0.9098798    0.2826087
## Max.cor.Y.rcv.3X3.glmnet           0.5962443    0.9098798    0.2826087
## Max.cor.Y.rcv.3X5.glmnet           0.5962443    0.9098798    0.2826087
## Max.cor.Y.rcv.5X1.glmnet           0.5962443    0.9098798    0.2826087
## Max.cor.Y.rcv.5X3.glmnet           0.5962443    0.9098798    0.2826087
## Max.cor.Y.rcv.5X5.glmnet           0.5962443    0.9098798    0.2826087
## Max.cor.Y.rcv.1X1.cp.0.rpart       0.6174697    0.9218959    0.3130435
## Max.cor.Y.rpart                    0.5870523    0.9045394    0.2695652
## Interact.High.cor.Y.glmnet         0.5962443    0.9098798    0.2826087
## Low.cor.X.glmnet                   0.5962443    0.9098798    0.2826087
## All.X.glmnet                       0.5965780    0.9105474    0.2826087
##                              max.AUCROCR.OOB opt.prob.threshold.OOB
## MFO.myMFO_classfr                  0.5000000                    0.1
## Random.myrandom_classfr            0.4857956                    0.1
## Max.cor.Y.rcv.1X1.glmnet           0.8116126                    0.1
## Max.cor.Y.rcv.3X1.glmnet           0.8067975                    0.1
## Max.cor.Y.rcv.3X3.glmnet           0.8067975                    0.1
## Max.cor.Y.rcv.3X5.glmnet           0.8067975                    0.1
## Max.cor.Y.rcv.5X1.glmnet           0.8114863                    0.1
## Max.cor.Y.rcv.5X3.glmnet           0.8114863                    0.1
## Max.cor.Y.rcv.5X5.glmnet           0.8114863                    0.1
## Max.cor.Y.rcv.1X1.cp.0.rpart       0.7773858                    0.1
## Max.cor.Y.rpart                    0.5892132                    0.6
## Interact.High.cor.Y.glmnet         0.8067975                    0.1
## Low.cor.X.glmnet                   0.8067975                    0.1
## All.X.glmnet                       0.8061009                    0.1
##                              max.f.score.OOB max.Accuracy.OOB
## MFO.myMFO_classfr                  0.2349336        0.1331019
## Random.myrandom_classfr            0.2349336        0.1331019
## Max.cor.Y.rcv.1X1.glmnet           0.4405405        0.7604167
## Max.cor.Y.rcv.3X1.glmnet           0.4375839        0.7575231
## Max.cor.Y.rcv.3X3.glmnet           0.4375839        0.7575231
## Max.cor.Y.rcv.3X5.glmnet           0.4375839        0.7575231
## Max.cor.Y.rcv.5X1.glmnet           0.4609375        0.7604167
## Max.cor.Y.rcv.5X3.glmnet           0.4609375        0.7604167
## Max.cor.Y.rcv.5X5.glmnet           0.4609375        0.7604167
## Max.cor.Y.rcv.1X1.cp.0.rpart       0.4207493        0.7673611
## Max.cor.Y.rpart                    0.2850575        0.8200231
## Interact.High.cor.Y.glmnet         0.4375839        0.7575231
## Low.cor.X.glmnet                   0.4375839        0.7575231
## All.X.glmnet                       0.4340426        0.7690972
##                              max.AccuracyLower.OOB max.AccuracyUpper.OOB
## MFO.myMFO_classfr                        0.1174298             0.1500310
## Random.myrandom_classfr                  0.1174298             0.1500310
## Max.cor.Y.rcv.1X1.glmnet                 0.7395703             0.7803749
## Max.cor.Y.rcv.3X1.glmnet                 0.7365992             0.7775689
## Max.cor.Y.rcv.3X3.glmnet                 0.7365992             0.7775689
## Max.cor.Y.rcv.3X5.glmnet                 0.7365992             0.7775689
## Max.cor.Y.rcv.5X1.glmnet                 0.7395703             0.7803749
## Max.cor.Y.rcv.5X3.glmnet                 0.7395703             0.7803749
## Max.cor.Y.rcv.5X5.glmnet                 0.7395703             0.7803749
## Max.cor.Y.rcv.1X1.cp.0.rpart             0.7467059             0.7871043
## Max.cor.Y.rpart                          0.8010821             0.8378705
## Interact.High.cor.Y.glmnet               0.7365992             0.7775689
## Low.cor.X.glmnet                         0.7365992             0.7775689
## All.X.glmnet                             0.7484910             0.7887855
##                              max.Kappa.OOB max.AccuracySD.fit
## MFO.myMFO_classfr                0.0000000                 NA
## Random.myrandom_classfr          0.0000000                 NA
## Max.cor.Y.rcv.1X1.glmnet         0.3148374                 NA
## Max.cor.Y.rcv.3X1.glmnet         0.3107477        0.007015493
## Max.cor.Y.rcv.3X3.glmnet         0.3107477        0.005178375
## Max.cor.Y.rcv.3X5.glmnet         0.3107477        0.005396525
## Max.cor.Y.rcv.5X1.glmnet         0.3373693        0.008837283
## Max.cor.Y.rcv.5X3.glmnet         0.3373693        0.006138477
## Max.cor.Y.rcv.5X5.glmnet         0.3373693        0.006213800
## Max.cor.Y.rcv.1X1.cp.0.rpart     0.2953321                 NA
## Max.cor.Y.rpart                  0.1825002        0.005069520
## Interact.High.cor.Y.glmnet       0.3107477        0.005178375
## Low.cor.X.glmnet                 0.3107477        0.005178375
## All.X.glmnet                     0.3103487        0.004892081
##                              max.KappaSD.fit
## MFO.myMFO_classfr                         NA
## Random.myrandom_classfr                   NA
## Max.cor.Y.rcv.1X1.glmnet                  NA
## Max.cor.Y.rcv.3X1.glmnet          0.02403706
## Max.cor.Y.rcv.3X3.glmnet          0.01754365
## Max.cor.Y.rcv.3X5.glmnet          0.01835474
## Max.cor.Y.rcv.5X1.glmnet          0.03133449
## Max.cor.Y.rcv.5X3.glmnet          0.02161286
## Max.cor.Y.rcv.5X5.glmnet          0.02210061
## Max.cor.Y.rcv.1X1.cp.0.rpart              NA
## Max.cor.Y.rpart                   0.01910910
## Interact.High.cor.Y.glmnet        0.01754365
## Low.cor.X.glmnet                  0.01754365
## All.X.glmnet                      0.01659075
```

```r
rm(ret_lst)
fit.models_1_chunk_df <- 
    myadd_chunk(fit.models_1_chunk_df, "fit.models_1_end", major.inc = FALSE,
                label.minor = "teardown")
```

```
##                  label step_major step_minor label_minor     bgn    end
## 4 fit.models_1_preProc          1          3     preProc 161.253 161.32
## 5     fit.models_1_end          1          4    teardown 161.320     NA
##   elapsed
## 4   0.067
## 5      NA
```

```r
glb_chunks_df <- myadd_chunk(glb_chunks_df, "fit.models", major.inc = FALSE)
```

```
##         label step_major step_minor label_minor     bgn     end elapsed
## 11 fit.models          6          1           1 139.591 161.329  21.738
## 12 fit.models          6          2           2 161.330      NA      NA
```


```r
fit.models_2_chunk_df <- 
    myadd_chunk(NULL, "fit.models_2_bgn", label.minor = "setup")
```

```
##              label step_major step_minor label_minor     bgn end elapsed
## 1 fit.models_2_bgn          1          0       setup 163.004  NA      NA
```

```r
#stop(here"); glb_to_sav(); all.equal(glb_models_df, sav_models_df)
# if (!is.null(glbMdlMetricSummaryFn)) {
#     stats_df <- glb_models_df[, "id", FALSE]
# 
#     stats_mdl_df <- data.frame()
#     for (mdl_id in stats_df$id) {
#         stats_mdl_df <- rbind(stats_mdl_df, 
#             mypredict_mdl(glb_models_lst[[mdl_id]], glbObsFit, glb_rsp_var, 
#                           glb_rsp_var_out, mdl_id, "fit",
#         						glbMdlMetricSummaryFn, glbMdlMetricSummary, 
#         						glbMdlMetricMaximize, ret_type="stats"))
#     }
#     stats_df <- merge(stats_df, stats_mdl_df, all.x=TRUE)
#     
#     stats_mdl_df <- data.frame()
#     for (mdl_id in stats_df$id) {
#         stats_mdl_df <- rbind(stats_mdl_df, 
#             mypredict_mdl(glb_models_lst[[mdl_id]], glbObsOOB, glb_rsp_var, 
#                           glb_rsp_var_out, mdl_id, "OOB",
#             					glbMdlMetricSummaryFn, glbMdlMetricSummary, 
#         						glbMdlMetricMaximize, ret_type="stats"))
#     }
#     stats_df <- merge(stats_df, stats_mdl_df, all.x=TRUE)
#     
#     print("Merging following data into glb_models_df:")
#     print(stats_mrg_df <- stats_df[, c(1, grep(glbMdlMetricSummary, names(stats_df)))])
#     print(tmp_models_df <- orderBy(~id, glb_models_df[, c("id",
#                                     grep(glbMdlMetricSummary, names(stats_df), value=TRUE))]))
# 
#     tmp2_models_df <- glb_models_df[, c("id", setdiff(names(glb_models_df),
#                                     grep(glbMdlMetricSummary, names(stats_df), value=TRUE)))]
#     tmp3_models_df <- merge(tmp2_models_df, stats_mrg_df, all.x=TRUE, sort=FALSE)
#     print(tmp3_models_df)
#     print(names(tmp3_models_df))
#     print(glb_models_df <- subset(tmp3_models_df, select=-id.1))
# }

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
##                                                        id
## MFO.myMFO_classfr                       MFO.myMFO_classfr
## Random.myrandom_classfr           Random.myrandom_classfr
## Max.cor.Y.rcv.1X1.glmnet         Max.cor.Y.rcv.1X1.glmnet
## Max.cor.Y.rcv.3X1.glmnet         Max.cor.Y.rcv.3X1.glmnet
## Max.cor.Y.rcv.3X3.glmnet         Max.cor.Y.rcv.3X3.glmnet
## Max.cor.Y.rcv.3X5.glmnet         Max.cor.Y.rcv.3X5.glmnet
## Max.cor.Y.rcv.5X1.glmnet         Max.cor.Y.rcv.5X1.glmnet
## Max.cor.Y.rcv.5X3.glmnet         Max.cor.Y.rcv.5X3.glmnet
## Max.cor.Y.rcv.5X5.glmnet         Max.cor.Y.rcv.5X5.glmnet
## Max.cor.Y.rcv.1X1.cp.0.rpart Max.cor.Y.rcv.1X1.cp.0.rpart
## Max.cor.Y.rpart                           Max.cor.Y.rpart
## Interact.High.cor.Y.glmnet     Interact.High.cor.Y.glmnet
## Low.cor.X.glmnet                         Low.cor.X.glmnet
## All.X.glmnet                                 All.X.glmnet
##                                                                                               feats
## MFO.myMFO_classfr                                                                            .rnorm
## Random.myrandom_classfr                                                                      .rnorm
## Max.cor.Y.rcv.1X1.glmnet                                           WordCount.root2,NDSSName.my.fctr
## Max.cor.Y.rcv.3X1.glmnet                                           WordCount.root2,NDSSName.my.fctr
## Max.cor.Y.rcv.3X3.glmnet                                           WordCount.root2,NDSSName.my.fctr
## Max.cor.Y.rcv.3X5.glmnet                                           WordCount.root2,NDSSName.my.fctr
## Max.cor.Y.rcv.5X1.glmnet                                           WordCount.root2,NDSSName.my.fctr
## Max.cor.Y.rcv.5X3.glmnet                                           WordCount.root2,NDSSName.my.fctr
## Max.cor.Y.rcv.5X5.glmnet                                           WordCount.root2,NDSSName.my.fctr
## Max.cor.Y.rcv.1X1.cp.0.rpart                                       WordCount.root2,NDSSName.my.fctr
## Max.cor.Y.rpart                                                    WordCount.root2,NDSSName.my.fctr
## Interact.High.cor.Y.glmnet         WordCount.root2,NDSSName.my.fctr,WordCount.root2:WordCount.root2
## Low.cor.X.glmnet                             WordCount.root2,NDSSName.my.fctr,.rnorm,WordCount.nexp
## All.X.glmnet                 WordCount.root2,WordCount.log1p,NDSSName.my.fctr,.rnorm,WordCount.nexp
##                              max.nTuningRuns max.AUCpROC.fit max.Sens.fit
## MFO.myMFO_classfr                          0       0.5000000    1.0000000
## Random.myrandom_classfr                    0       0.4990604    0.8312611
## Max.cor.Y.rcv.1X1.glmnet                   0       0.8790544    0.9632073
## Max.cor.Y.rcv.3X1.glmnet                  25       0.8767919    0.9644760
## Max.cor.Y.rcv.3X3.glmnet                  25       0.8767919    0.9644760
## Max.cor.Y.rcv.3X5.glmnet                  25       0.8767919    0.9644760
## Max.cor.Y.rcv.5X1.glmnet                  25       0.8784031    0.9642223
## Max.cor.Y.rcv.5X3.glmnet                  25       0.8784031    0.9642223
## Max.cor.Y.rcv.5X5.glmnet                  25       0.8784031    0.9642223
## Max.cor.Y.rcv.1X1.cp.0.rpart               0       0.8821543    0.9705658
## Max.cor.Y.rpart                            5       0.8709432    0.9632073
## Interact.High.cor.Y.glmnet                25       0.8767919    0.9644760
## Low.cor.X.glmnet                          25       0.8767919    0.9644760
## All.X.glmnet                              25       0.8771175    0.9639685
##                              max.Spec.fit max.AUCROCR.fit
## MFO.myMFO_classfr               0.0000000       0.5000000
## Random.myrandom_classfr         0.1668598       0.4972757
## Max.cor.Y.rcv.1X1.glmnet        0.7949015       0.9608594
## Max.cor.Y.rcv.3X1.glmnet        0.7891078       0.9582555
## Max.cor.Y.rcv.3X3.glmnet        0.7891078       0.9582555
## Max.cor.Y.rcv.3X5.glmnet        0.7891078       0.9582555
## Max.cor.Y.rcv.5X1.glmnet        0.7925840       0.9607052
## Max.cor.Y.rcv.5X3.glmnet        0.7925840       0.9607052
## Max.cor.Y.rcv.5X5.glmnet        0.7925840       0.9607052
## Max.cor.Y.rcv.1X1.cp.0.rpart    0.7937428       0.9504198
## Max.cor.Y.rpart                 0.7786790       0.8746354
## Interact.High.cor.Y.glmnet      0.7891078       0.9582555
## Low.cor.X.glmnet                0.7891078       0.9582555
## All.X.glmnet                    0.7902665       0.9589072
##                              opt.prob.threshold.fit max.f.score.fit
## MFO.myMFO_classfr                               0.1       0.3045703
## Random.myrandom_classfr                         0.1       0.3045703
## Max.cor.Y.rcv.1X1.glmnet                        0.5       0.8099174
## Max.cor.Y.rcv.3X1.glmnet                        0.4       0.8099174
## Max.cor.Y.rcv.3X3.glmnet                        0.4       0.8099174
## Max.cor.Y.rcv.3X5.glmnet                        0.4       0.8099174
## Max.cor.Y.rcv.5X1.glmnet                        0.5       0.8104265
## Max.cor.Y.rcv.5X3.glmnet                        0.5       0.8104265
## Max.cor.Y.rcv.5X5.glmnet                        0.5       0.8104265
## Max.cor.Y.rcv.1X1.cp.0.rpart                    0.4       0.8235294
## Max.cor.Y.rpart                                 0.6       0.8000000
## Interact.High.cor.Y.glmnet                      0.4       0.8099174
## Low.cor.X.glmnet                                0.4       0.8099174
## All.X.glmnet                                    0.4       0.8091873
##                              max.Accuracy.fit max.Kappa.fit
## MFO.myMFO_classfr                   0.1796420     0.0000000
## Random.myrandom_classfr             0.1796420     0.0000000
## Max.cor.Y.rcv.1X1.glmnet            0.9329725     0.7692476
## Max.cor.Y.rcv.3X1.glmnet            0.9335973     0.7691678
## Max.cor.Y.rcv.3X3.glmnet            0.9333193     0.7690803
## Max.cor.Y.rcv.3X5.glmnet            0.9332218     0.7686375
## Max.cor.Y.rcv.5X1.glmnet            0.9331818     0.7689055
## Max.cor.Y.rcv.5X3.glmnet            0.9333905     0.7698577
## Max.cor.Y.rcv.5X5.glmnet            0.9331816     0.7691429
## Max.cor.Y.rcv.1X1.cp.0.rpart        0.9381765     0.7860827
## Max.cor.Y.rpart                     0.9296422     0.7515134
## Interact.High.cor.Y.glmnet          0.9333193     0.7690803
## Low.cor.X.glmnet                    0.9333193     0.7690803
## All.X.glmnet                        0.9326952     0.7668828
##                              max.AUCpROC.OOB max.Sens.OOB max.Spec.OOB
## MFO.myMFO_classfr                  0.5000000    1.0000000    0.0000000
## Random.myrandom_classfr            0.5125675    0.8077437    0.2173913
## Max.cor.Y.rcv.1X1.glmnet           0.5962443    0.9098798    0.2826087
## Max.cor.Y.rcv.3X1.glmnet           0.5962443    0.9098798    0.2826087
## Max.cor.Y.rcv.3X3.glmnet           0.5962443    0.9098798    0.2826087
## Max.cor.Y.rcv.3X5.glmnet           0.5962443    0.9098798    0.2826087
## Max.cor.Y.rcv.5X1.glmnet           0.5962443    0.9098798    0.2826087
## Max.cor.Y.rcv.5X3.glmnet           0.5962443    0.9098798    0.2826087
## Max.cor.Y.rcv.5X5.glmnet           0.5962443    0.9098798    0.2826087
## Max.cor.Y.rcv.1X1.cp.0.rpart       0.6174697    0.9218959    0.3130435
## Max.cor.Y.rpart                    0.5870523    0.9045394    0.2695652
## Interact.High.cor.Y.glmnet         0.5962443    0.9098798    0.2826087
## Low.cor.X.glmnet                   0.5962443    0.9098798    0.2826087
## All.X.glmnet                       0.5965780    0.9105474    0.2826087
##                              max.AUCROCR.OOB opt.prob.threshold.OOB
## MFO.myMFO_classfr                  0.5000000                    0.1
## Random.myrandom_classfr            0.4857956                    0.1
## Max.cor.Y.rcv.1X1.glmnet           0.8116126                    0.1
## Max.cor.Y.rcv.3X1.glmnet           0.8067975                    0.1
## Max.cor.Y.rcv.3X3.glmnet           0.8067975                    0.1
## Max.cor.Y.rcv.3X5.glmnet           0.8067975                    0.1
## Max.cor.Y.rcv.5X1.glmnet           0.8114863                    0.1
## Max.cor.Y.rcv.5X3.glmnet           0.8114863                    0.1
## Max.cor.Y.rcv.5X5.glmnet           0.8114863                    0.1
## Max.cor.Y.rcv.1X1.cp.0.rpart       0.7773858                    0.1
## Max.cor.Y.rpart                    0.5892132                    0.6
## Interact.High.cor.Y.glmnet         0.8067975                    0.1
## Low.cor.X.glmnet                   0.8067975                    0.1
## All.X.glmnet                       0.8061009                    0.1
##                              max.f.score.OOB max.Accuracy.OOB
## MFO.myMFO_classfr                  0.2349336        0.1331019
## Random.myrandom_classfr            0.2349336        0.1331019
## Max.cor.Y.rcv.1X1.glmnet           0.4405405        0.7604167
## Max.cor.Y.rcv.3X1.glmnet           0.4375839        0.7575231
## Max.cor.Y.rcv.3X3.glmnet           0.4375839        0.7575231
## Max.cor.Y.rcv.3X5.glmnet           0.4375839        0.7575231
## Max.cor.Y.rcv.5X1.glmnet           0.4609375        0.7604167
## Max.cor.Y.rcv.5X3.glmnet           0.4609375        0.7604167
## Max.cor.Y.rcv.5X5.glmnet           0.4609375        0.7604167
## Max.cor.Y.rcv.1X1.cp.0.rpart       0.4207493        0.7673611
## Max.cor.Y.rpart                    0.2850575        0.8200231
## Interact.High.cor.Y.glmnet         0.4375839        0.7575231
## Low.cor.X.glmnet                   0.4375839        0.7575231
## All.X.glmnet                       0.4340426        0.7690972
##                              max.Kappa.OOB inv.elapsedtime.everything
## MFO.myMFO_classfr                0.0000000                  3.7174721
## Random.myrandom_classfr          0.0000000                  3.7174721
## Max.cor.Y.rcv.1X1.glmnet         0.3148374                  1.0288066
## Max.cor.Y.rcv.3X1.glmnet         0.3107477                  0.3538570
## Max.cor.Y.rcv.3X3.glmnet         0.3107477                  0.2075981
## Max.cor.Y.rcv.3X5.glmnet         0.3107477                  0.1557875
## Max.cor.Y.rcv.5X1.glmnet         0.3373693                  0.3102699
## Max.cor.Y.rcv.5X3.glmnet         0.3373693                  0.1628134
## Max.cor.Y.rcv.5X5.glmnet         0.3373693                  0.1122586
## Max.cor.Y.rcv.1X1.cp.0.rpart     0.2953321                  1.1587486
## Max.cor.Y.rpart                  0.1825002                  0.3469813
## Interact.High.cor.Y.glmnet       0.3107477                  0.2275313
## Low.cor.X.glmnet                 0.3107477                  0.2198285
## All.X.glmnet                     0.3103487                  0.2075981
##                              inv.elapsedtime.final
## MFO.myMFO_classfr                       333.333333
## Random.myrandom_classfr                1000.000000
## Max.cor.Y.rcv.1X1.glmnet                  3.663004
## Max.cor.Y.rcv.3X1.glmnet                  3.703704
## Max.cor.Y.rcv.3X3.glmnet                  3.703704
## Max.cor.Y.rcv.3X5.glmnet                  3.717472
## Max.cor.Y.rcv.5X1.glmnet                  3.745318
## Max.cor.Y.rcv.5X3.glmnet                  3.745318
## Max.cor.Y.rcv.5X5.glmnet                  3.731343
## Max.cor.Y.rcv.1X1.cp.0.rpart             14.285714
## Max.cor.Y.rpart                          14.084507
## Interact.High.cor.Y.glmnet                3.703704
## Low.cor.X.glmnet                          3.412969
## All.X.glmnet                              3.401361
```

```r
print(myplot_radar(radar_inp_df=plt_models_df))
```

```
## Warning in RColorBrewer::brewer.pal(n, pal): n too large, allowed maximum for palette Set1 is 9
## Returning the palette you asked for with that many colors
```

```
## Warning: The shape palette can deal with a maximum of 6 discrete values
## because more than 6 becomes difficult to discriminate; you have
## 14. Consider specifying shapes manually if you must have them.
```

```
## Warning: Removed 160 rows containing missing values (geom_point).
```

```
## Warning in RColorBrewer::brewer.pal(n, pal): n too large, allowed maximum for palette Set1 is 9
## Returning the palette you asked for with that many colors
```

```
## Warning: The shape palette can deal with a maximum of 6 discrete values
## because more than 6 becomes difficult to discriminate; you have
## 14. Consider specifying shapes manually if you must have them.
```

![](NYTBlogs3_base_files/figure-html/fit.models_2-1.png) 

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
## Warning: Removed 4 rows containing missing values (geom_errorbar).
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
## Warning: Removed 4 rows containing missing values (geom_errorbar).
```

![](NYTBlogs3_base_files/figure-html/fit.models_2-2.png) 

```r
dsp_models_cols <- c("id", 
                    glbMdlMetricsEval[glbMdlMetricsEval %in% names(glb_models_df)],
                    grep("opt.", names(glb_models_df), fixed = TRUE, value = TRUE)) 
# if (glb_is_classification && glb_is_binomial) 
#     dsp_models_cols <- c(dsp_models_cols, "opt.prob.threshold.OOB")
print(dsp_models_df <- orderBy(get_model_sel_frmla(), glb_models_df)[, dsp_models_cols])
```

```
##                                                        id max.Accuracy.OOB
## Max.cor.Y.rpart                           Max.cor.Y.rpart        0.8200231
## All.X.glmnet                                 All.X.glmnet        0.7690972
## Max.cor.Y.rcv.1X1.cp.0.rpart Max.cor.Y.rcv.1X1.cp.0.rpart        0.7673611
## Max.cor.Y.rcv.1X1.glmnet         Max.cor.Y.rcv.1X1.glmnet        0.7604167
## Max.cor.Y.rcv.5X3.glmnet         Max.cor.Y.rcv.5X3.glmnet        0.7604167
## Max.cor.Y.rcv.5X1.glmnet         Max.cor.Y.rcv.5X1.glmnet        0.7604167
## Max.cor.Y.rcv.5X5.glmnet         Max.cor.Y.rcv.5X5.glmnet        0.7604167
## Max.cor.Y.rcv.3X1.glmnet         Max.cor.Y.rcv.3X1.glmnet        0.7575231
## Max.cor.Y.rcv.3X3.glmnet         Max.cor.Y.rcv.3X3.glmnet        0.7575231
## Interact.High.cor.Y.glmnet     Interact.High.cor.Y.glmnet        0.7575231
## Low.cor.X.glmnet                         Low.cor.X.glmnet        0.7575231
## Max.cor.Y.rcv.3X5.glmnet         Max.cor.Y.rcv.3X5.glmnet        0.7575231
## MFO.myMFO_classfr                       MFO.myMFO_classfr        0.1331019
## Random.myrandom_classfr           Random.myrandom_classfr        0.1331019
##                              max.AUCROCR.OOB max.AUCpROC.OOB
## Max.cor.Y.rpart                    0.5892132       0.5870523
## All.X.glmnet                       0.8061009       0.5965780
## Max.cor.Y.rcv.1X1.cp.0.rpart       0.7773858       0.6174697
## Max.cor.Y.rcv.1X1.glmnet           0.8116126       0.5962443
## Max.cor.Y.rcv.5X3.glmnet           0.8114863       0.5962443
## Max.cor.Y.rcv.5X1.glmnet           0.8114863       0.5962443
## Max.cor.Y.rcv.5X5.glmnet           0.8114863       0.5962443
## Max.cor.Y.rcv.3X1.glmnet           0.8067975       0.5962443
## Max.cor.Y.rcv.3X3.glmnet           0.8067975       0.5962443
## Interact.High.cor.Y.glmnet         0.8067975       0.5962443
## Low.cor.X.glmnet                   0.8067975       0.5962443
## Max.cor.Y.rcv.3X5.glmnet           0.8067975       0.5962443
## MFO.myMFO_classfr                  0.5000000       0.5000000
## Random.myrandom_classfr            0.4857956       0.5125675
##                              max.Accuracy.fit opt.prob.threshold.fit
## Max.cor.Y.rpart                     0.9296422                    0.6
## All.X.glmnet                        0.9326952                    0.4
## Max.cor.Y.rcv.1X1.cp.0.rpart        0.9381765                    0.4
## Max.cor.Y.rcv.1X1.glmnet            0.9329725                    0.5
## Max.cor.Y.rcv.5X3.glmnet            0.9333905                    0.5
## Max.cor.Y.rcv.5X1.glmnet            0.9331818                    0.5
## Max.cor.Y.rcv.5X5.glmnet            0.9331816                    0.5
## Max.cor.Y.rcv.3X1.glmnet            0.9335973                    0.4
## Max.cor.Y.rcv.3X3.glmnet            0.9333193                    0.4
## Interact.High.cor.Y.glmnet          0.9333193                    0.4
## Low.cor.X.glmnet                    0.9333193                    0.4
## Max.cor.Y.rcv.3X5.glmnet            0.9332218                    0.4
## MFO.myMFO_classfr                   0.1796420                    0.1
## Random.myrandom_classfr             0.1796420                    0.1
##                              opt.prob.threshold.OOB
## Max.cor.Y.rpart                                 0.6
## All.X.glmnet                                    0.1
## Max.cor.Y.rcv.1X1.cp.0.rpart                    0.1
## Max.cor.Y.rcv.1X1.glmnet                        0.1
## Max.cor.Y.rcv.5X3.glmnet                        0.1
## Max.cor.Y.rcv.5X1.glmnet                        0.1
## Max.cor.Y.rcv.5X5.glmnet                        0.1
## Max.cor.Y.rcv.3X1.glmnet                        0.1
## Max.cor.Y.rcv.3X3.glmnet                        0.1
## Interact.High.cor.Y.glmnet                      0.1
## Low.cor.X.glmnet                                0.1
## Max.cor.Y.rcv.3X5.glmnet                        0.1
## MFO.myMFO_classfr                               0.1
## Random.myrandom_classfr                         0.1
```

```r
print(myplot_radar(radar_inp_df = dsp_models_df))
```

```
## Warning in RColorBrewer::brewer.pal(n, pal): n too large, allowed maximum for palette Set1 is 9
## Returning the palette you asked for with that many colors
```

```
## Warning: The shape palette can deal with a maximum of 6 discrete values
## because more than 6 becomes difficult to discriminate; you have
## 14. Consider specifying shapes manually if you must have them.
```

```
## Warning: Removed 56 rows containing missing values (geom_point).
```

```
## Warning in RColorBrewer::brewer.pal(n, pal): n too large, allowed maximum for palette Set1 is 9
## Returning the palette you asked for with that many colors
```

```
## Warning: The shape palette can deal with a maximum of 6 discrete values
## because more than 6 becomes difficult to discriminate; you have
## 14. Consider specifying shapes manually if you must have them.
```

![](NYTBlogs3_base_files/figure-html/fit.models_2-3.png) 

```r
print("Metrics used for model selection:"); print(get_model_sel_frmla())
```

```
## [1] "Metrics used for model selection:"
```

```
## ~-max.Accuracy.OOB - max.AUCROCR.OOB - max.AUCpROC.OOB - max.Accuracy.fit - 
##     opt.prob.threshold.OOB
## <environment: 0x7ff241329310>
```

```r
print(sprintf("Best model id: %s", dsp_models_df[1, "id"]))
```

```
## [1] "Best model id: Max.cor.Y.rpart"
```

```r
glb_get_predictions <- function(df, mdl_id, rsp_var_out, prob_threshold_def=NULL, verbose=FALSE) {
    mdl <- glb_models_lst[[mdl_id]]
    #rsp_var_out <- paste0(rsp_var_out, mdl_id)
    
    rsp_var_out <- paste0(glb_rsp_var, ".predict.")
    predct_var_name <- paste0(rsp_var_out, mdl_id)        
    predct_prob_var_name <- paste0(rsp_var_out, mdl_id, ".prob")    
    predct_accurate_var_name <- paste0(rsp_var_out, mdl_id, ".accurate")
    predct_error_var_name <- paste0(rsp_var_out, mdl_id, ".err")
    predct_erabs_var_name <- paste0(rsp_var_out, mdl_id, ".err.abs")

    if (glb_is_regression) {
        df[, predct_var_name] <- predict(mdl, newdata=df, type="raw")
        if (verbose) print(myplot_scatter(df, glb_rsp_var, predct_var_name) + 
                  facet_wrap(reformulate(glb_category_var), scales = "free") + 
                  stat_smooth(method="glm"))

        df[, predct_error_var_name] <- df[, predct_var_name] - df[, glb_rsp_var]
        if (verbose) print(myplot_scatter(df, predct_var_name, predct_error_var_name) + 
                  #facet_wrap(reformulate(glb_category_var), scales = "free") + 
                  stat_smooth(method="auto"))
        if (verbose) print(myplot_scatter(df, glb_rsp_var, predct_error_var_name) + 
                  #facet_wrap(reformulate(glb_category_var), scales = "free") + 
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
#                   facet_wrap(reformulate(glb_category_var), scales = "free") + 
#                   stat_smooth(method="glm"))

        df[, predct_error_var_name] <- df[, predct_var_name] != df[, glb_rsp_var]
#         if (verbose) print(myplot_scatter(df, predct_var_name, predct_error_var_name) + 
#                   #facet_wrap(reformulate(glb_category_var), scales = "free") + 
#                   stat_smooth(method="auto"))
#         if (verbose) print(myplot_scatter(df, glb_rsp_var, predct_error_var_name) + 
#                   #facet_wrap(reformulate(glb_category_var), scales = "free") + 
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

#stop(here"); glb_to_sav(); glbObsAll <- sav_allobs_df; glbObsTrn <- sav_trnobs_df; glbObsFit <- sav_fitobs_df; glbObsOOB <- sav_OOBobs_df; sav_models_df <- glb_models_df; glb_models_df <- sav_models_df; glb_featsimp_df <- sav_featsimp_df    

myget_category_stats <- function(obs_df, mdl_id, label) {
    require(dplyr)
    require(lazyeval)
    
    predct_var_name <- paste0(glb_rsp_var_out, mdl_id)        
    predct_error_var_name <- paste0(glb_rsp_var_out, mdl_id, ".err.abs")
    
    if (!predct_var_name %in% names(obs_df))
        obs_df <- glb_get_predictions(obs_df, mdl_id, glb_rsp_var_out)
    
    tmp_obs_df <- obs_df %>%
        dplyr::select_(glb_category_var, glb_rsp_var, predct_var_name, predct_error_var_name) 
    #dplyr::rename(startprice.log10.predict.RFE.X.glmnet.err=error_abs_OOB)
    names(tmp_obs_df)[length(names(tmp_obs_df))] <- paste0("err.abs.", label)
    
    ret_ctgry_df <- tmp_obs_df %>%
        dplyr::group_by_(glb_category_var) %>%
        dplyr::summarise_(#interp(~sum(abs(var)), var=as.name(glb_rsp_var)), 
            interp(~sum(var), var=as.name(paste0("err.abs.", label))), 
            interp(~mean(var), var=as.name(paste0("err.abs.", label))),
            interp(~n()))
    names(ret_ctgry_df) <- c(glb_category_var, 
                             #paste0(glb_rsp_var, ".abs.", label, ".sum"),
                             paste0("err.abs.", label, ".sum"),                             
                             paste0("err.abs.", label, ".mean"), 
                             paste0(".n.", label))
    ret_ctgry_df <- dplyr::ungroup(ret_ctgry_df)
    #colSums(ret_ctgry_df[, -grep(glb_category_var, names(ret_ctgry_df))])
    
    return(ret_ctgry_df)    
}
#print(colSums((ctgry_df <- myget_category_stats(obs_df=glbObsFit, mdl_id="", label="fit"))[, -grep(glb_category_var, names(ctgry_df))]))

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
        glbObsFit <- glb_get_predictions(df = glbObsFit, mdl_id,
                                             glb_rsp_var_out)
        glbObsOOB <- glb_get_predictions(df = glbObsOOB, mdl_id,
                                             glb_rsp_var_out)
    }
    
#mdl_id_pfx <- "Ensemble.RFE"; mdlId <- paste0(mdl_id_pfx, ".glmnet")
#glb_mdl_ensemble <- gsub(glb_rsp_var_out, "", grep("RFE\\.X\\.(?!Interact)", row.names(glb_featsimp_df), perl = TRUE, value = TRUE), fixed = TRUE)
#varImp(glb_models_lst[[mdlId]])
    
#cor_df <- data.frame(cor=cor(glbObsFit[, glb_rsp_var], glbObsFit[, paste(glb_rsp_var_out, glb_mdl_ensemble)], use="pairwise.complete.obs"))
#glbObsFit <- glb_get_predictions(df=glbObsFit, "Ensemble.glmnet", glb_rsp_var_out);print(colSums((ctgry_df <- myget_category_stats(obs_df=glbObsFit, mdl_id="Ensemble.glmnet", label="fit"))[, -grep(glb_category_var, names(ctgry_df))]))
    
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

    indep_vars <- paste(glb_rsp_var_out, glb_mdl_ensemble, sep = "")
    if (glb_is_classification)
        indep_vars <- paste(indep_vars, ".prob", sep = "")
    # Some models in glb_mdl_ensemble might not be fitted e.g. RFE.X.Interact
    indep_vars <- intersect(indep_vars, names(glbObsFit))
    
#     indep_vars <- grep(glb_rsp_var_out, names(glbObsFit), fixed=TRUE, value=TRUE)
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
## [1] "User specified selection: All.X.glmnet"
```

```r
myprint_mdl(glb_sel_mdl <- glb_models_lst[[glb_sel_mdl_id]])
```

![](NYTBlogs3_base_files/figure-html/fit.models_2-4.png) 

```
##             Length Class      Mode     
## a0            91   -none-     numeric  
## beta        2184   dgCMatrix  S4       
## df            91   -none-     numeric  
## dim            2   -none-     numeric  
## lambda        91   -none-     numeric  
## dev.ratio     91   -none-     numeric  
## nulldev        1   -none-     numeric  
## npasses        1   -none-     numeric  
## jerr           1   -none-     numeric  
## offset         1   -none-     logical  
## classnames     2   -none-     character
## call           5   -none-     call     
## nobs           1   -none-     numeric  
## lambdaOpt      1   -none-     numeric  
## xNames        24   -none-     character
## problemType    1   -none-     character
## tuneValue      2   data.frame list     
## obsLevels      2   -none-     character
## [1] "min lambda > lambdaOpt:"
##                                        (Intercept) 
##                                        -5.52332827 
##                       NDSSName.my.fctr#Multimedia# 
##                                        -0.16474095 
##             NDSSName.my.fctr#Opinion#RoomForDebate 
##                                        -2.18132553 
##           NDSSName.my.fctr#Opinion#ThePublicEditor 
##                                         3.90787932 
##                    NDSSName.my.fctr#U.S.#Education 
##                                        -0.90764328 
## NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness 
##                                        -0.24502067 
##         NDSSName.my.fctrBusiness#Crosswords/Games# 
##                                         4.67941504 
##               NDSSName.my.fctrBusiness#Technology# 
##                                         0.83419892 
##                     NDSSName.my.fctrForeign#World# 
##                                        -0.08805297 
##          NDSSName.my.fctrForeign#World#AsiaPacific 
##                                        -0.83659624 
##                      NDSSName.my.fctrOpEd#Opinion# 
##                                         4.79265681 
##                    NDSSName.my.fctrScience#Health# 
##                                         3.67959505 
##                    NDSSName.my.fctrStyles##Fashion 
##                                        -0.64461547 
##                       NDSSName.my.fctrStyles#U.S.# 
##                                         3.44343706 
##                           NDSSName.my.fctrTStyle## 
##                                        -0.97236910 
##                            NDSSName.my.fctrmyOther 
##                                        -0.23762498 
##                                    WordCount.log1p 
##                                         0.18619880 
##                                    WordCount.root2 
##                                         0.06691104 
## [1] "max lambda < lambdaOpt:"
##                                        (Intercept) 
##                                        -5.59881022 
##                       NDSSName.my.fctr#Multimedia# 
##                                        -0.24609678 
##             NDSSName.my.fctr#Opinion#RoomForDebate 
##                                        -2.30818954 
##           NDSSName.my.fctr#Opinion#ThePublicEditor 
##                                         3.94273791 
##                    NDSSName.my.fctr#U.S.#Education 
##                                        -1.01466873 
## NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness 
##                                        -0.30479819 
##         NDSSName.my.fctrBusiness#Crosswords/Games# 
##                                         4.69839363 
##               NDSSName.my.fctrBusiness#Technology# 
##                                         0.84455935 
##                     NDSSName.my.fctrForeign#World# 
##                                        -0.18204983 
##          NDSSName.my.fctrForeign#World#AsiaPacific 
##                                        -0.91223472 
##                      NDSSName.my.fctrOpEd#Opinion# 
##                                         4.80649138 
##                    NDSSName.my.fctrScience#Health# 
##                                         3.68561250 
##                    NDSSName.my.fctrStyles##Fashion 
##                                        -0.74721366 
##                       NDSSName.my.fctrStyles#U.S.# 
##                                         3.44963924 
##                           NDSSName.my.fctrTStyle## 
##                                        -1.04099014 
##                     NDSSName.my.fctrTravel#Travel# 
##                                        -0.03203556 
##                            NDSSName.my.fctrmyOther 
##                                        -0.35136209 
##                                    WordCount.log1p 
##                                         0.19889748 
##                                    WordCount.root2 
##                                         0.06700663
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
## [1] "All.X.glmnet fit prediction diagnostics:"
```

```r
glbObsFit <- glb_get_predictions(df=glbObsFit, mdl_id=glb_sel_mdl_id, 
                                     rsp_var_out=glb_rsp_var_out)
print(sprintf("%s OOB prediction diagnostics:", glb_sel_mdl_id))
```

```
## [1] "All.X.glmnet OOB prediction diagnostics:"
```

```r
glbObsOOB <- glb_get_predictions(df = glbObsOOB, mdl_id = glb_sel_mdl_id, 
                                     rsp_var_out = glb_rsp_var_out)

glb_featsimp_df <- 
    myget_feats_importance(mdl=glb_sel_mdl, featsimp_df=NULL)
glb_featsimp_df[, paste0(glb_sel_mdl_id, ".importance")] <- glb_featsimp_df$importance
#mdl_id <-"RFE.X.glmnet"; glb_featsimp_df <- myget_feats_importance(glb_models_lst[[mdl_id]], glb_featsimp_df); glb_featsimp_df[, paste0(mdl_id, ".importance")] <- glb_featsimp_df$importance; print(glb_featsimp_df)
#print(head(sbst_featsimp_df <- subset(glb_featsimp_df, is.na(RFE.X.glmnet.importance) | (abs(RFE.X.YeoJohnson.glmnet.importance - RFE.X.glmnet.importance) > 0.0001), select=-importance)))
#print(orderBy(~ -cor.y.abs, subset(glb_feats_df, id %in% c(row.names(sbst_featsimp_df), "startprice.dcm1.is9", "D.weight.post.stop.sum"))))
print(glb_featsimp_df)
```

```
##                                                    importance
## NDSSName.my.fctrOpEd#Opinion#                       100.00000
## NDSSName.my.fctrBusiness#Crosswords/Games#           98.38173
## NDSSName.my.fctr#Opinion#ThePublicEditor             87.34197
## NDSSName.my.fctrScience#Health#                      84.05064
## NDSSName.my.fctrStyles#U.S.#                         80.66804
## NDSSName.my.fctrBusiness#Technology#                 43.29623
## WordCount.log1p                                      34.01597
## WordCount.root2                                      32.29795
## .rnorm                                               31.33944
## NDSSName.my.fctrBusiness#BusinessDay#Dealbook        31.33944
## NDSSName.my.fctrCulture##                            31.33944
## NDSSName.my.fctrCulture#Arts#                        31.33944
## NDSSName.my.fctrMetro#N.Y./Region#                   31.33944
## WordCount.nexp                                       31.33944
## NDSSName.my.fctrTravel#Travel#                       31.31570
## NDSSName.my.fctrForeign#World#                       30.00852
## NDSSName.my.fctr#Multimedia#                         28.91940
## NDSSName.my.fctrmyOther                              27.85141
## NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness   27.78546
## NDSSName.my.fctrStyles##Fashion                      22.02991
## NDSSName.my.fctrForeign#World#AsiaPacific            19.29994
## NDSSName.my.fctr#U.S.#Education                      18.25900
## NDSSName.my.fctrTStyle##                             17.36032
## NDSSName.my.fctr#Opinion#RoomForDebate                0.00000
##                                                    All.X.glmnet.importance
## NDSSName.my.fctrOpEd#Opinion#                                    100.00000
## NDSSName.my.fctrBusiness#Crosswords/Games#                        98.38173
## NDSSName.my.fctr#Opinion#ThePublicEditor                          87.34197
## NDSSName.my.fctrScience#Health#                                   84.05064
## NDSSName.my.fctrStyles#U.S.#                                      80.66804
## NDSSName.my.fctrBusiness#Technology#                              43.29623
## WordCount.log1p                                                   34.01597
## WordCount.root2                                                   32.29795
## .rnorm                                                            31.33944
## NDSSName.my.fctrBusiness#BusinessDay#Dealbook                     31.33944
## NDSSName.my.fctrCulture##                                         31.33944
## NDSSName.my.fctrCulture#Arts#                                     31.33944
## NDSSName.my.fctrMetro#N.Y./Region#                                31.33944
## WordCount.nexp                                                    31.33944
## NDSSName.my.fctrTravel#Travel#                                    31.31570
## NDSSName.my.fctrForeign#World#                                    30.00852
## NDSSName.my.fctr#Multimedia#                                      28.91940
## NDSSName.my.fctrmyOther                                           27.85141
## NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness                27.78546
## NDSSName.my.fctrStyles##Fashion                                   22.02991
## NDSSName.my.fctrForeign#World#AsiaPacific                         19.29994
## NDSSName.my.fctr#U.S.#Education                                   18.25900
## NDSSName.my.fctrTStyle##                                          17.36032
## NDSSName.my.fctr#Opinion#RoomForDebate                             0.00000
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
        featsimp_df <- orderBy(~ -importance.max, 
            summaryBy(importance ~ feat + feat.interact, data=featsimp_df,
                      FUN=max))    
        #rex_str=":(.*)"; txt_vctr=tail(featsimp_df$feat); ret_lst <- regexec(rex_str, txt_vctr); ret_lst <- regmatches(txt_vctr, ret_lst); ret_vctr <- sapply(1:length(ret_lst), function(pos_ix) ifelse(length(ret_lst[[pos_ix]]) > 0, ret_lst[[pos_ix]], "")); print(ret_vctr <- ret_vctr[ret_vctr != ""])    
        
        featsimp_df <- subset(featsimp_df, !is.na(importance.max))
        if (nrow(featsimp_df) > 5) {
            warning("Limiting important feature scatter plots to 5 out of ",
                    nrow(featsimp_df))
            featsimp_df <- head(featsimp_df, 5)
        }
        
    #     if (!all(is.na(featsimp_df$feat.interact)))
    #         stop("not implemented yet")
        rsp_var_out <- paste0(glb_rsp_var_out, mdl_id)
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
    glb_analytics_diag_plots(obs_df=glbObsOOB, mdl_id=glb_sel_mdl_id)                  
```

![](NYTBlogs3_base_files/figure-html/fit.models_2-5.png) ![](NYTBlogs3_base_files/figure-html/fit.models_2-6.png) ![](NYTBlogs3_base_files/figure-html/fit.models_2-7.png) ![](NYTBlogs3_base_files/figure-html/fit.models_2-8.png) ![](NYTBlogs3_base_files/figure-html/fit.models_2-9.png) 

```
## [1] "Min/Max Boundaries: "
##      UniqueID Popular.fctr Popular.fctr.predict.All.X.glmnet.prob
## 1618     1618            N                            0.003961226
## 6251     6251            Y                            0.933618142
## 6435     6435            N                            0.110713384
##      Popular.fctr.predict.All.X.glmnet
## 1618                                 N
## 6251                                 Y
## 6435                                 Y
##      Popular.fctr.predict.All.X.glmnet.err
## 1618                                 FALSE
## 6251                                 FALSE
## 6435                                  TRUE
##      Popular.fctr.predict.All.X.glmnet.err.abs
## 1618                               0.003961226
## 6251                               0.066381858
## 6435                               0.110713384
##      Popular.fctr.predict.All.X.glmnet.accurate
## 1618                                       TRUE
## 6251                                       TRUE
## 6435                                      FALSE
##      Popular.fctr.predict.All.X.glmnet.error .label
## 1618                              0.00000000   1618
## 6251                              0.00000000   6251
## 6435                              0.01071338   6435
## [1] "Inaccurate: "
##      UniqueID Popular.fctr Popular.fctr.predict.All.X.glmnet.prob
## 4020     4020            Y                             0.01153216
## 4352     4352            Y                             0.01557858
## 4775     4775            Y                             0.01694357
## 6354     6354            Y                             0.01942833
## 4745     4745            Y                             0.02014592
## 172       172            Y                             0.02244134
##      Popular.fctr.predict.All.X.glmnet
## 4020                                 N
## 4352                                 N
## 4775                                 N
## 6354                                 N
## 4745                                 N
## 172                                  N
##      Popular.fctr.predict.All.X.glmnet.err
## 4020                                  TRUE
## 4352                                  TRUE
## 4775                                  TRUE
## 6354                                  TRUE
## 4745                                  TRUE
## 172                                   TRUE
##      Popular.fctr.predict.All.X.glmnet.err.abs
## 4020                                 0.9884678
## 4352                                 0.9844214
## 4775                                 0.9830564
## 6354                                 0.9805717
## 4745                                 0.9798541
## 172                                  0.9775587
##      Popular.fctr.predict.All.X.glmnet.accurate
## 4020                                      FALSE
## 4352                                      FALSE
## 4775                                      FALSE
## 6354                                      FALSE
## 4745                                      FALSE
## 172                                       FALSE
##      Popular.fctr.predict.All.X.glmnet.error
## 4020                             -0.08846784
## 4352                             -0.08442142
## 4775                             -0.08305643
## 6354                             -0.08057167
## 4745                             -0.07985408
## 172                              -0.07755866
##      UniqueID Popular.fctr Popular.fctr.predict.All.X.glmnet.prob
## 2545     2545            N                              0.1014865
## 4917     4917            N                              0.1072512
## 303       303            N                              0.1511238
## 5087     5087            N                              0.5422170
## 4962     4962            N                              0.6157458
## 6511     6511            N                              0.8246298
##      Popular.fctr.predict.All.X.glmnet
## 2545                                 Y
## 4917                                 Y
## 303                                  Y
## 5087                                 Y
## 4962                                 Y
## 6511                                 Y
##      Popular.fctr.predict.All.X.glmnet.err
## 2545                                  TRUE
## 4917                                  TRUE
## 303                                   TRUE
## 5087                                  TRUE
## 4962                                  TRUE
## 6511                                  TRUE
##      Popular.fctr.predict.All.X.glmnet.err.abs
## 2545                                 0.1014865
## 4917                                 0.1072512
## 303                                  0.1511238
## 5087                                 0.5422170
## 4962                                 0.6157458
## 6511                                 0.8246298
##      Popular.fctr.predict.All.X.glmnet.accurate
## 2545                                      FALSE
## 4917                                      FALSE
## 303                                       FALSE
## 5087                                      FALSE
## 4962                                      FALSE
## 6511                                      FALSE
##      Popular.fctr.predict.All.X.glmnet.error
## 2545                             0.001486507
## 4917                             0.007251173
## 303                              0.051123801
## 5087                             0.442216987
## 4962                             0.515745834
## 6511                             0.724629843
##      UniqueID Popular.fctr Popular.fctr.predict.All.X.glmnet.prob
## 770       770            N                              0.9698459
## 6276     6276            N                              0.9721252
## 221       221            N                              0.9729676
## 3590     3590            N                              0.9730589
## 472       472            N                              0.9740792
## 2995     2995            N                              0.9744222
##      Popular.fctr.predict.All.X.glmnet
## 770                                  Y
## 6276                                 Y
## 221                                  Y
## 3590                                 Y
## 472                                  Y
## 2995                                 Y
##      Popular.fctr.predict.All.X.glmnet.err
## 770                                   TRUE
## 6276                                  TRUE
## 221                                   TRUE
## 3590                                  TRUE
## 472                                   TRUE
## 2995                                  TRUE
##      Popular.fctr.predict.All.X.glmnet.err.abs
## 770                                  0.9698459
## 6276                                 0.9721252
## 221                                  0.9729676
## 3590                                 0.9730589
## 472                                  0.9740792
## 2995                                 0.9744222
##      Popular.fctr.predict.All.X.glmnet.accurate
## 770                                       FALSE
## 6276                                      FALSE
## 221                                       FALSE
## 3590                                      FALSE
## 472                                       FALSE
## 2995                                      FALSE
##      Popular.fctr.predict.All.X.glmnet.error
## 770                                0.8698459
## 6276                               0.8721252
## 221                                0.8729676
## 3590                               0.8730589
## 472                                0.8740792
## 2995                               0.8744222
```

![](NYTBlogs3_base_files/figure-html/fit.models_2-10.png) 

```r
if (!is.null(glb_category_var)) {
    glb_ctgry_df <- merge(glb_ctgry_df, 
            myget_category_stats(obs_df = glbObsFit, mdl_id = glb_sel_mdl_id, label = "fit"),
                          by = glb_category_var, all = TRUE)
    row.names(glb_ctgry_df) <- glb_ctgry_df[, glb_category_var]
    glb_ctgry_df <- merge(glb_ctgry_df, 
                myget_category_stats(obs_df=glbObsOOB, mdl_id=glb_sel_mdl_id, label="OOB"),
                          #by=glb_category_var, all=TRUE) glb_ctgry-df already contains .n.OOB ?
                          all=TRUE)
    row.names(glb_ctgry_df) <- glb_ctgry_df[, glb_category_var]
    if (any(grepl("OOB", glbMdlMetricsEval)))
        print(orderBy(~-err.abs.OOB.mean, glb_ctgry_df)) else
            print(orderBy(~-err.abs.fit.mean, glb_ctgry_df))
    print(colSums(glb_ctgry_df[, -grep(glb_category_var, names(glb_ctgry_df))]))
}
```

```
##                                                      NDSSName.my.fctr
## OpEd#Opinion#                                           OpEd#Opinion#
## Styles#U.S.#                                             Styles#U.S.#
## Science#Health#                                       Science#Health#
## Business#Crosswords/Games#                 Business#Crosswords/Games#
## #Opinion#ThePublicEditor                     #Opinion#ThePublicEditor
## Business#Technology#                             Business#Technology#
## Business#BusinessDay#Dealbook           Business#BusinessDay#Dealbook
## ##                                                                 ##
## Metro#N.Y./Region#                                 Metro#N.Y./Region#
## Culture#Arts#                                           Culture#Arts#
## Business#BusinessDay#SmallBusiness Business#BusinessDay#SmallBusiness
## Styles##Fashion                                       Styles##Fashion
## #Opinion#RoomForDebate                         #Opinion#RoomForDebate
## #Multimedia#                                             #Multimedia#
## Travel#Travel#                                         Travel#Travel#
## TStyle##                                                     TStyle##
## Foreign#World#AsiaPacific                   Foreign#World#AsiaPacific
## myOther                                                       myOther
## Culture##                                                   Culture##
## #U.S.#Education                                       #U.S.#Education
## Foreign#World#                                         Foreign#World#
##                                    .n.OOB .n.Fit .n.Tst .freqRatio.Fit
## OpEd#Opinion#                          89    437    164    0.090965862
## Styles#U.S.#                           50    127     61    0.026436303
## Science#Health#                        48    148     57    0.030807660
## Business#Crosswords/Games#             18    105     42    0.021856786
## #Opinion#ThePublicEditor                4     16     10    0.003330558
## Business#Technology#                  126    213    114    0.044338052
## Business#BusinessDay#Dealbook         323    629    304    0.130932556
## ##                                    371    913    342    0.190049958
## Metro#N.Y./Region#                     70    128     67    0.026644463
## Culture#Arts#                         185    490    174    0.101998335
## Business#BusinessDay#SmallBusiness     40    100     41    0.020815987
## Styles##Fashion                        15    104     15    0.021648626
## #Opinion#RoomForDebate                 20     42     20    0.008742714
## #Multimedia#                           49     92     52    0.019150708
## Travel#Travel#                         34     83     35    0.017277269
## TStyle##                              101    623    105    0.129683597
## Foreign#World#AsiaPacific              53    150     56    0.031223980
## myOther                                 5     33      5    0.006869276
## Culture##                               1     NA     70             NA
## #U.S.#Education                        82    243     89    0.050582848
## Foreign#World#                         44    128     47    0.026644463
##                                    .freqRatio.OOB .freqRatio.Tst
## OpEd#Opinion#                        0.0515046296    0.087700535
## Styles#U.S.#                         0.0289351852    0.032620321
## Science#Health#                      0.0277777778    0.030481283
## Business#Crosswords/Games#           0.0104166667    0.022459893
## #Opinion#ThePublicEditor             0.0023148148    0.005347594
## Business#Technology#                 0.0729166667    0.060962567
## Business#BusinessDay#Dealbook        0.1869212963    0.162566845
## ##                                   0.2146990741    0.182887701
## Metro#N.Y./Region#                   0.0405092593    0.035828877
## Culture#Arts#                        0.1070601852    0.093048128
## Business#BusinessDay#SmallBusiness   0.0231481481    0.021925134
## Styles##Fashion                      0.0086805556    0.008021390
## #Opinion#RoomForDebate               0.0115740741    0.010695187
## #Multimedia#                         0.0283564815    0.027807487
## Travel#Travel#                       0.0196759259    0.018716578
## TStyle##                             0.0584490741    0.056149733
## Foreign#World#AsiaPacific            0.0306712963    0.029946524
## myOther                              0.0028935185    0.002673797
## Culture##                            0.0005787037    0.037433155
## #U.S.#Education                      0.0474537037    0.047593583
## Foreign#World#                       0.0254629630    0.025133690
##                                    err.abs.fit.sum err.abs.fit.mean .n.fit
## OpEd#Opinion#                           101.962286       0.23332331    437
## Styles#U.S.#                             51.707342       0.40714443    127
## Science#Health#                          48.666469       0.32882749    148
## Business#Crosswords/Games#               19.346059       0.18424818    105
## #Opinion#ThePublicEditor                  4.044658       0.25279116     16
## Business#Technology#                     46.446343       0.21805795    213
## Business#BusinessDay#Dealbook            76.983498       0.12239030    629
## ##                                       91.828632       0.10057901    913
## Metro#N.Y./Region#                       14.444993       0.11285150    128
## Culture#Arts#                            41.599488       0.08489691    490
## Business#BusinessDay#SmallBusiness        8.607239       0.08607239    100
## Styles##Fashion                           3.033172       0.02916512    104
## #Opinion#RoomForDebate                    1.938007       0.04614303     42
## #Multimedia#                              4.739265       0.05151375     92
## Travel#Travel#                            2.569132       0.03095340     83
## TStyle##                                 14.962143       0.02401628    623
## Foreign#World#AsiaPacific                 7.489991       0.04993328    150
## myOther                                   1.720179       0.05212665     33
## Culture##                                       NA               NA     NA
## #U.S.#Education                           4.569704       0.01880537    243
## Foreign#World#                            3.357443       0.02623002    128
##                                    err.abs.OOB.sum err.abs.OOB.mean
## OpEd#Opinion#                          52.11236621       0.58553220
## Styles#U.S.#                           27.39861838       0.54797237
## Science#Health#                        25.12604617       0.52345930
## Business#Crosswords/Games#              9.01410297       0.50078350
## #Opinion#ThePublicEditor                1.97530353       0.49382588
## Business#Technology#                   28.55232120       0.22660572
## Business#BusinessDay#Dealbook          57.75141149       0.17879694
## ##                                     63.39672286       0.17088065
## Metro#N.Y./Region#                     11.12676976       0.15895385
## Culture#Arts#                          28.98681899       0.15668551
## Business#BusinessDay#SmallBusiness      3.84496527       0.09612413
## Styles##Fashion                         1.42711318       0.09514088
## #Opinion#RoomForDebate                  1.78910542       0.08945527
## #Multimedia#                            2.85858104       0.05833839
## Travel#Travel#                          1.97769164       0.05816740
## TStyle##                                5.15834308       0.05107270
## Foreign#World#AsiaPacific               2.69695832       0.05088601
## myOther                                 0.23180569       0.04636114
## Culture##                               0.02849447       0.02849447
## #U.S.#Education                         2.22191005       0.02709646
## Foreign#World#                          1.13916144       0.02589003
##           .n.OOB           .n.Fit           .n.Tst   .freqRatio.Fit 
##      1728.000000               NA      1870.000000               NA 
##   .freqRatio.OOB   .freqRatio.Tst  err.abs.fit.sum err.abs.fit.mean 
##         1.000000         1.000000               NA               NA 
##           .n.fit  err.abs.OOB.sum err.abs.OOB.mean 
##               NA       328.814611         4.170523
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
## 1 fit.models_2_bgn          1          0    teardown 174.883  NA      NA
```

```r
glb_chunks_df <- myadd_chunk(glb_chunks_df, "fit.models", major.inc=FALSE)
```

```
##         label step_major step_minor label_minor     bgn     end elapsed
## 12 fit.models          6          2           2 161.330 174.893  13.563
## 13 fit.models          6          3           3 174.894      NA      NA
```


```r
# if (sum(is.na(glbObsAll$D.P.http)) > 0)
#         stop("fit.models_3: Why is this happening ?")

#stop(here"); glb_to_sav()
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

![](NYTBlogs3_base_files/figure-html/fit.models_3-1.png) 

```r
glb_chunks_df <- myadd_chunk(glb_chunks_df, "fit.data.training", major.inc=TRUE)
```

```
##                label step_major step_minor label_minor     bgn     end
## 13        fit.models          6          3           3 174.894 181.343
## 14 fit.data.training          7          0           0 181.344      NA
##    elapsed
## 13    6.45
## 14      NA
```

## Step `7.0: fit data training`

```r
#load(paste0(glb_inp_pfx, "dsk.RData"))

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
        mdlimp_df <- subset(myget_feats_importance(glb_sel_mdl), importance > 5)
        # Fit selected models on glbObsTrn
        for (mdl_id in gsub(".prob", "", 
                    gsub(glb_rsp_var_out, "", row.names(mdlimp_df), fixed = TRUE),
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
                        type = glb_model_type, tune.df = glb_tune_models_df,
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
                                                rsp_var_out = glb_rsp_var_out,
                                                prob_threshold_def = 
                    subset(glb_models_df, id == mdl_id)$opt.prob.threshold.OOB)
            glbObsNew <- glb_get_predictions(df = glbObsNew,
                                                mdl_id = tail(glb_models_df$id, 1), 
                                                rsp_var_out = glb_rsp_var_out,
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
        mdlimp_df <- subset(myget_feats_importance(glb_sel_mdl), importance > 5)
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
        
    # Discontinuing use of tune_finmdl_df; 
    #   since final model needs to be cved on glbObsTrn
#     tune_finmdl_df <- NULL
#     if (nrow(glb_sel_mdl$bestTune) > 0) {
#         for (param in names(glb_sel_mdl$bestTune)) {
#             #print(sprintf("param: %s", param))
#             if (glb_sel_mdl$bestTune[1, param] != "none")
#                 tune_finmdl_df <- rbind(tune_finmdl_df, 
#                     data.frame(parameter=param, 
#                                min=glb_sel_mdl$bestTune[1, param], 
#                                max=glb_sel_mdl$bestTune[1, param], 
#                                by=1)) # by val does not matter
#         }
#     } 
    
    # Sync with parameters in mydsutils.R
#stop(here"); glb_to_sav(); glb_models_lst <- sav_models_lst; glb_models_df <- sav_models_df
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
    method_vctr <- 
        unique(c("glm", tail(unlist(strsplit(glb_sel_mdl_id, "[.]")), 1)))
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
trainControl.allowParallel = FALSE, 
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
## +(rfe) fit Fold1.Rep1 size: 23 
## -(rfe) fit Fold1.Rep1 size: 23 
## +(rfe) imp Fold1.Rep1 
## -(rfe) imp Fold1.Rep1 
## +(rfe) fit Fold1.Rep1 size: 16 
## -(rfe) fit Fold1.Rep1 size: 16 
## +(rfe) fit Fold1.Rep1 size:  8 
## -(rfe) fit Fold1.Rep1 size:  8 
## +(rfe) fit Fold1.Rep1 size:  4 
## -(rfe) fit Fold1.Rep1 size:  4 
## +(rfe) fit Fold1.Rep1 size:  2 
## -(rfe) fit Fold1.Rep1 size:  2 
## +(rfe) fit Fold2.Rep1 size: 23 
## -(rfe) fit Fold2.Rep1 size: 23 
## +(rfe) imp Fold2.Rep1 
## -(rfe) imp Fold2.Rep1 
## +(rfe) fit Fold2.Rep1 size: 16 
## -(rfe) fit Fold2.Rep1 size: 16 
## +(rfe) fit Fold2.Rep1 size:  8 
## -(rfe) fit Fold2.Rep1 size:  8 
## +(rfe) fit Fold2.Rep1 size:  4 
## -(rfe) fit Fold2.Rep1 size:  4 
## +(rfe) fit Fold2.Rep1 size:  2 
## -(rfe) fit Fold2.Rep1 size:  2 
## +(rfe) fit Fold3.Rep1 size: 23 
## -(rfe) fit Fold3.Rep1 size: 23 
## +(rfe) imp Fold3.Rep1 
## -(rfe) imp Fold3.Rep1 
## +(rfe) fit Fold3.Rep1 size: 16 
## -(rfe) fit Fold3.Rep1 size: 16 
## +(rfe) fit Fold3.Rep1 size:  8 
## -(rfe) fit Fold3.Rep1 size:  8 
## +(rfe) fit Fold3.Rep1 size:  4 
## -(rfe) fit Fold3.Rep1 size:  4 
## +(rfe) fit Fold3.Rep1 size:  2 
## -(rfe) fit Fold3.Rep1 size:  2 
## +(rfe) fit Fold1.Rep2 size: 23 
## -(rfe) fit Fold1.Rep2 size: 23 
## +(rfe) imp Fold1.Rep2 
## -(rfe) imp Fold1.Rep2 
## +(rfe) fit Fold1.Rep2 size: 16 
## -(rfe) fit Fold1.Rep2 size: 16 
## +(rfe) fit Fold1.Rep2 size:  8 
## -(rfe) fit Fold1.Rep2 size:  8 
## +(rfe) fit Fold1.Rep2 size:  4 
## -(rfe) fit Fold1.Rep2 size:  4 
## +(rfe) fit Fold1.Rep2 size:  2 
## -(rfe) fit Fold1.Rep2 size:  2 
## +(rfe) fit Fold2.Rep2 size: 23 
## -(rfe) fit Fold2.Rep2 size: 23 
## +(rfe) imp Fold2.Rep2 
## -(rfe) imp Fold2.Rep2 
## +(rfe) fit Fold2.Rep2 size: 16 
## -(rfe) fit Fold2.Rep2 size: 16 
## +(rfe) fit Fold2.Rep2 size:  8 
## -(rfe) fit Fold2.Rep2 size:  8 
## +(rfe) fit Fold2.Rep2 size:  4 
## -(rfe) fit Fold2.Rep2 size:  4 
## +(rfe) fit Fold2.Rep2 size:  2 
## -(rfe) fit Fold2.Rep2 size:  2 
## +(rfe) fit Fold3.Rep2 size: 23 
## -(rfe) fit Fold3.Rep2 size: 23 
## +(rfe) imp Fold3.Rep2 
## -(rfe) imp Fold3.Rep2 
## +(rfe) fit Fold3.Rep2 size: 16 
## -(rfe) fit Fold3.Rep2 size: 16 
## +(rfe) fit Fold3.Rep2 size:  8 
## -(rfe) fit Fold3.Rep2 size:  8 
## +(rfe) fit Fold3.Rep2 size:  4 
## -(rfe) fit Fold3.Rep2 size:  4 
## +(rfe) fit Fold3.Rep2 size:  2 
## -(rfe) fit Fold3.Rep2 size:  2 
## +(rfe) fit Fold1.Rep3 size: 23 
## -(rfe) fit Fold1.Rep3 size: 23 
## +(rfe) imp Fold1.Rep3 
## -(rfe) imp Fold1.Rep3 
## +(rfe) fit Fold1.Rep3 size: 16 
## -(rfe) fit Fold1.Rep3 size: 16 
## +(rfe) fit Fold1.Rep3 size:  8 
## -(rfe) fit Fold1.Rep3 size:  8 
## +(rfe) fit Fold1.Rep3 size:  4 
## -(rfe) fit Fold1.Rep3 size:  4 
## +(rfe) fit Fold1.Rep3 size:  2 
## -(rfe) fit Fold1.Rep3 size:  2 
## +(rfe) fit Fold2.Rep3 size: 23 
## -(rfe) fit Fold2.Rep3 size: 23 
## +(rfe) imp Fold2.Rep3 
## -(rfe) imp Fold2.Rep3 
## +(rfe) fit Fold2.Rep3 size: 16 
## -(rfe) fit Fold2.Rep3 size: 16 
## +(rfe) fit Fold2.Rep3 size:  8 
## -(rfe) fit Fold2.Rep3 size:  8 
## +(rfe) fit Fold2.Rep3 size:  4 
## -(rfe) fit Fold2.Rep3 size:  4 
## +(rfe) fit Fold2.Rep3 size:  2 
## -(rfe) fit Fold2.Rep3 size:  2 
## +(rfe) fit Fold3.Rep3 size: 23 
## -(rfe) fit Fold3.Rep3 size: 23 
## +(rfe) imp Fold3.Rep3 
## -(rfe) imp Fold3.Rep3 
## +(rfe) fit Fold3.Rep3 size: 16 
## -(rfe) fit Fold3.Rep3 size: 16 
## +(rfe) fit Fold3.Rep3 size:  8 
## -(rfe) fit Fold3.Rep3 size:  8 
## +(rfe) fit Fold3.Rep3 size:  4 
## -(rfe) fit Fold3.Rep3 size:  4 
## +(rfe) fit Fold3.Rep3 size:  2 
## -(rfe) fit Fold3.Rep3 size:  2 
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
##          8   0.9018 0.63845   0.006472 0.02305         
##         16   0.9016 0.63798   0.006364 0.02240         
##         23   0.9033 0.64712   0.006639 0.02346        *
## 
## The top 5 variables (out of 23):
##    WordCount.log1p, WordCount.root2, WordCount.nexp, NDSSName.my.fctrOpEd#Opinion#, NDSSName.my.fctrScience#Health#
## 
##  [1] "WordCount.log1p"                                   
##  [2] "WordCount.root2"                                   
##  [3] "WordCount.nexp"                                    
##  [4] "NDSSName.my.fctrOpEd#Opinion#"                     
##  [5] "NDSSName.my.fctrScience#Health#"                   
##  [6] "NDSSName.my.fctrBusiness#Crosswords/Games#"        
##  [7] "NDSSName.my.fctrStyles#U.S.#"                      
##  [8] ".rnorm"                                            
##  [9] "NDSSName.my.fctrBusiness#Technology#"              
## [10] "NDSSName.my.fctrmyOther"                           
## [11] "NDSSName.my.fctr#Opinion#RoomForDebate"            
## [12] "NDSSName.my.fctrMetro#N.Y./Region#"                
## [13] "NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness"
## [14] "NDSSName.my.fctrTravel#Travel#"                    
## [15] "NDSSName.my.fctrStyles##Fashion"                   
## [16] "NDSSName.my.fctr#Multimedia#"                      
## [17] "NDSSName.my.fctrForeign#World#"                    
## [18] "NDSSName.my.fctrForeign#World#AsiaPacific"         
## [19] "NDSSName.my.fctr#U.S.#Education"                   
## [20] "NDSSName.my.fctrCulture#Arts#"                     
## [21] "NDSSName.my.fctrBusiness#BusinessDay#Dealbook"     
## [22] "NDSSName.my.fctr##"                                
## [23] "NDSSName.my.fctrTStyle##"
```

```
## [1] "fitting model: Final.glm"
## [1] "    indep_vars: WordCount.root2,WordCount.log1p,NDSSName.my.fctr,.rnorm,WordCount.nexp"
## + Fold1.Rep1: parameter=none
```

```
## Warning: glm.fit: algorithm did not converge
```

```
## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
```

```
## - Fold1.Rep1: parameter=none 
## + Fold2.Rep1: parameter=none
```

```
## Warning in predict.lm(object, newdata, se.fit, scale = 1, type =
## ifelse(type == : prediction from a rank-deficient fit may be misleading
```

```
## - Fold2.Rep1: parameter=none 
## + Fold3.Rep1: parameter=none 
## - Fold3.Rep1: parameter=none 
## + Fold1.Rep2: parameter=none 
## - Fold1.Rep2: parameter=none 
## + Fold2.Rep2: parameter=none 
## - Fold2.Rep2: parameter=none 
## + Fold3.Rep2: parameter=none
```

```
## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
```

```
## Warning in predict.lm(object, newdata, se.fit, scale = 1, type =
## ifelse(type == : prediction from a rank-deficient fit may be misleading
```

```
## - Fold3.Rep2: parameter=none 
## + Fold1.Rep3: parameter=none
```

```
## Warning in predict.lm(object, newdata, se.fit, scale = 1, type =
## ifelse(type == : prediction from a rank-deficient fit may be misleading
```

```
## - Fold1.Rep3: parameter=none 
## + Fold2.Rep3: parameter=none 
## - Fold2.Rep3: parameter=none 
## + Fold3.Rep3: parameter=none
```

```
## Warning: glm.fit: algorithm did not converge
```

```
## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
```

```
## - Fold3.Rep3: parameter=none 
## Aggregating results
## Fitting final model on full training set
```

```
## Warning: not plotting observations with leverage one:
##   4746
```

![](NYTBlogs3_base_files/figure-html/fit.data.training_0-1.png) ![](NYTBlogs3_base_files/figure-html/fit.data.training_0-2.png) ![](NYTBlogs3_base_files/figure-html/fit.data.training_0-3.png) 

```
## Warning: not plotting observations with leverage one:
##   4746
```

![](NYTBlogs3_base_files/figure-html/fit.data.training_0-4.png) ![](NYTBlogs3_base_files/figure-html/fit.data.training_0-5.png) 

```
## 
## Call:
## NULL
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -2.6763  -0.3646  -0.2216  -0.0040   3.4318  
## 
## Coefficients:
##                                                        Estimate Std. Error
## (Intercept)                                          -6.375e+00  7.803e-01
## .rnorm                                               -8.854e-03  4.645e-02
## `NDSSName.my.fctr#Multimedia#`                       -1.972e+00  7.249e-01
## `NDSSName.my.fctr#Opinion#RoomForDebate`             -4.099e+00  1.029e+00
## `NDSSName.my.fctr#Opinion#ThePublicEditor`            3.212e+00  5.811e-01
## `NDSSName.my.fctr#U.S.#Education`                    -1.623e+01  3.335e+02
## `NDSSName.my.fctrBusiness#BusinessDay#Dealbook`      -4.406e-01  1.572e-01
## `NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness` -1.444e+00  4.697e-01
## `NDSSName.my.fctrBusiness#Crosswords/Games#`          3.786e+00  2.732e-01
## `NDSSName.my.fctrBusiness#Technology#`                4.470e-01  1.896e-01
## `NDSSName.my.fctrCulture##`                          -1.549e+01  6.523e+03
## `NDSSName.my.fctrCulture#Arts#`                      -5.450e-02  1.876e-01
## `NDSSName.my.fctrForeign#World#`                     -1.550e+01  4.776e+02
## `NDSSName.my.fctrForeign#World#AsiaPacific`          -2.497e+00  5.941e-01
## `NDSSName.my.fctrMetro#N.Y./Region#`                 -1.490e-01  2.852e-01
## `NDSSName.my.fctrOpEd#Opinion#`                       3.742e+00  1.616e-01
## `NDSSName.my.fctrScience#Health#`                     2.670e+00  1.905e-01
## `NDSSName.my.fctrStyles##Fashion`                    -2.572e+00  1.011e+00
## `NDSSName.my.fctrStyles#U.S.#`                        2.358e+00  1.952e-01
## `NDSSName.my.fctrTStyle##`                           -2.135e+00  3.738e-01
## `NDSSName.my.fctrTravel#Travel#`                     -1.816e+00  1.012e+00
## NDSSName.my.fctrmyOther                              -1.650e+01  9.845e+02
## WordCount.log1p                                       4.480e-01  1.827e-01
## WordCount.nexp                                       -4.904e+00  2.796e+01
## WordCount.root2                                       6.029e-02  1.631e-02
##                                                      z value Pr(>|z|)    
## (Intercept)                                           -8.170 3.09e-16 ***
## .rnorm                                                -0.191 0.848831    
## `NDSSName.my.fctr#Multimedia#`                        -2.720 0.006521 ** 
## `NDSSName.my.fctr#Opinion#RoomForDebate`              -3.985 6.76e-05 ***
## `NDSSName.my.fctr#Opinion#ThePublicEditor`             5.528 3.24e-08 ***
## `NDSSName.my.fctr#U.S.#Education`                     -0.049 0.961192    
## `NDSSName.my.fctrBusiness#BusinessDay#Dealbook`       -2.803 0.005068 ** 
## `NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness`  -3.075 0.002103 ** 
## `NDSSName.my.fctrBusiness#Crosswords/Games#`          13.862  < 2e-16 ***
## `NDSSName.my.fctrBusiness#Technology#`                 2.358 0.018384 *  
## `NDSSName.my.fctrCulture##`                           -0.002 0.998106    
## `NDSSName.my.fctrCulture#Arts#`                       -0.290 0.771462    
## `NDSSName.my.fctrForeign#World#`                      -0.032 0.974115    
## `NDSSName.my.fctrForeign#World#AsiaPacific`           -4.204 2.63e-05 ***
## `NDSSName.my.fctrMetro#N.Y./Region#`                  -0.522 0.601428    
## `NDSSName.my.fctrOpEd#Opinion#`                       23.151  < 2e-16 ***
## `NDSSName.my.fctrScience#Health#`                     14.015  < 2e-16 ***
## `NDSSName.my.fctrStyles##Fashion`                     -2.544 0.010969 *  
## `NDSSName.my.fctrStyles#U.S.#`                        12.084  < 2e-16 ***
## `NDSSName.my.fctrTStyle##`                            -5.710 1.13e-08 ***
## `NDSSName.my.fctrTravel#Travel#`                      -1.795 0.072576 .  
## NDSSName.my.fctrmyOther                               -0.017 0.986627    
## WordCount.log1p                                        2.452 0.014212 *  
## WordCount.nexp                                        -0.175 0.860772    
## WordCount.root2                                        3.698 0.000218 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 5900.1  on 6531  degrees of freedom
## Residual deviance: 3110.5  on 6507  degrees of freedom
## AIC: 3160.5
## 
## Number of Fisher Scoring iterations: 17
```

![](NYTBlogs3_base_files/figure-html/fit.data.training_0-6.png) ![](NYTBlogs3_base_files/figure-html/fit.data.training_0-7.png) 

```
##   Popular.fctr Popular.fctr.predict.Final.glm.N
## 1            N                             5172
## 2            Y                              363
##   Popular.fctr.predict.Final.glm.Y
## 1                              267
## 2                              730
##          Prediction
## Reference    N    Y
##         N 5172  267
##         Y  363  730
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   9.035517e-01   6.413003e-01   8.961350e-01   9.106059e-01   8.326699e-01 
## AccuracyPValue  McnemarPValue 
##   7.824873e-61   1.537762e-04 
##          id
## 1 Final.glm
##                                                                    feats
## 1 WordCount.root2,WordCount.log1p,NDSSName.my.fctr,.rnorm,WordCount.nexp
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1               1                      3.307                 0.261
##   max.AUCpROC.fit max.Sens.fit max.Spec.fit max.AUCROCR.fit
## 1       0.7989009    0.9610222    0.6367795       0.9333994
##   opt.prob.threshold.fit max.f.score.fit max.Accuracy.fit
## 1                    0.4       0.6985646        0.9057971
##   max.AccuracyLower.fit max.AccuracyUpper.fit max.Kappa.fit
## 1              0.896135             0.9106059     0.6371833
##   max.AccuracySD.fit max.KappaSD.fit
## 1        0.004986598      0.02341093
## [1] "fitting model: Final.glmnet"
## [1] "    indep_vars: WordCount.root2,WordCount.log1p,NDSSName.my.fctr,.rnorm,WordCount.nexp"
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
## Fitting alpha = 0.325, lambda = 0.0168 on full training set
```

![](NYTBlogs3_base_files/figure-html/fit.data.training_0-8.png) ![](NYTBlogs3_base_files/figure-html/fit.data.training_0-9.png) 

```
##             Length Class      Mode     
## a0            93   -none-     numeric  
## beta        2232   dgCMatrix  S4       
## df            93   -none-     numeric  
## dim            2   -none-     numeric  
## lambda        93   -none-     numeric  
## dev.ratio     93   -none-     numeric  
## nulldev        1   -none-     numeric  
## npasses        1   -none-     numeric  
## jerr           1   -none-     numeric  
## offset         1   -none-     logical  
## classnames     2   -none-     character
## call           5   -none-     call     
## nobs           1   -none-     numeric  
## lambdaOpt      1   -none-     numeric  
## xNames        24   -none-     character
## problemType    1   -none-     character
## tuneValue      2   data.frame list     
## obsLevels      2   -none-     character
## [1] "min lambda > lambdaOpt:"
##                                        (Intercept) 
##                                        -5.01773303 
##                       NDSSName.my.fctr#Multimedia# 
##                                        -0.31536157 
##             NDSSName.my.fctr#Opinion#RoomForDebate 
##                                        -1.51314492 
##           NDSSName.my.fctr#Opinion#ThePublicEditor 
##                                         2.67877297 
##                    NDSSName.my.fctr#U.S.#Education 
##                                        -0.82284989 
##      NDSSName.my.fctrBusiness#BusinessDay#Dealbook 
##                                        -0.02723099 
## NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness 
##                                        -0.37946894 
##         NDSSName.my.fctrBusiness#Crosswords/Games# 
##                                         3.34578807 
##               NDSSName.my.fctrBusiness#Technology# 
##                                         0.31843410 
##                     NDSSName.my.fctrForeign#World# 
##                                        -0.27776127 
##          NDSSName.my.fctrForeign#World#AsiaPacific 
##                                        -0.81472898 
##                      NDSSName.my.fctrOpEd#Opinion# 
##                                         3.35163474 
##                    NDSSName.my.fctrScience#Health# 
##                                         2.48264645 
##                    NDSSName.my.fctrStyles##Fashion 
##                                        -0.49241325 
##                       NDSSName.my.fctrStyles#U.S.# 
##                                         2.22613995 
##                           NDSSName.my.fctrTStyle## 
##                                        -0.72232021 
##                     NDSSName.my.fctrTravel#Travel# 
##                                        -0.07844427 
##                            NDSSName.my.fctrmyOther 
##                                        -0.22253379 
##                                    WordCount.log1p 
##                                         0.26612836 
##                                    WordCount.root2 
##                                         0.04648425 
## [1] "max lambda < lambdaOpt:"
##                                        (Intercept) 
##                                        -5.09261870 
##                       NDSSName.my.fctr#Multimedia# 
##                                        -0.37993993 
##             NDSSName.my.fctr#Opinion#RoomForDebate 
##                                        -1.62062746 
##           NDSSName.my.fctr#Opinion#ThePublicEditor 
##                                         2.72614928 
##                    NDSSName.my.fctr#U.S.#Education 
##                                        -0.89871970 
##      NDSSName.my.fctrBusiness#BusinessDay#Dealbook 
##                                        -0.04679393 
## NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness 
##                                        -0.43224773 
##         NDSSName.my.fctrBusiness#Crosswords/Games# 
##                                         3.38410449 
##               NDSSName.my.fctrBusiness#Technology# 
##                                         0.33911630 
##                     NDSSName.my.fctrForeign#World# 
##                                        -0.34241871 
##          NDSSName.my.fctrForeign#World#AsiaPacific 
##                                        -0.88251852 
##                      NDSSName.my.fctrOpEd#Opinion# 
##                                         3.38386039 
##                    NDSSName.my.fctrScience#Health# 
##                                         2.50486732 
##                    NDSSName.my.fctrStyles##Fashion 
##                                        -0.56144626 
##                       NDSSName.my.fctrStyles#U.S.# 
##                                         2.24594116 
##                           NDSSName.my.fctrTStyle## 
##                                        -0.77459454 
##                     NDSSName.my.fctrTravel#Travel# 
##                                        -0.13420086 
##                            NDSSName.my.fctrmyOther 
##                                        -0.31775172 
##                                    WordCount.log1p 
##                                         0.27449523 
##                                    WordCount.root2 
##                                         0.04746532
```

![](NYTBlogs3_base_files/figure-html/fit.data.training_0-10.png) ![](NYTBlogs3_base_files/figure-html/fit.data.training_0-11.png) 

```
##   Popular.fctr Popular.fctr.predict.Final.glmnet.N
## 1            N                                5141
## 2            Y                                 338
##   Popular.fctr.predict.Final.glmnet.Y
## 1                                 298
## 2                                 755
##          Prediction
## Reference    N    Y
##         N 5141  298
##         Y  338  755
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   9.026332e-01   6.454065e-01   8.951865e-01   9.097181e-01   8.326699e-01 
## AccuracyPValue  McnemarPValue 
##   3.401653e-59   1.219958e-01 
##             id
## 1 Final.glmnet
##                                                                    feats
## 1 WordCount.root2,WordCount.log1p,NDSSName.my.fctr,.rnorm,WordCount.nexp
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1              25                     19.047                 0.389
##   max.AUCpROC.fit max.Sens.fit max.Spec.fit max.AUCROCR.fit
## 1       0.7945167    0.9641478    0.6248856       0.9310545
##   opt.prob.threshold.fit max.f.score.fit max.Accuracy.fit
## 1                    0.3       0.7036347         0.907124
##   max.AccuracyLower.fit max.AccuracyUpper.fit max.Kappa.fit
## 1             0.8951865             0.9097181     0.6380888
##   max.AccuracySD.fit max.KappaSD.fit
## 1        0.004657963      0.02195699
```

```r
rm(ret_lst)
glb_chunks_df <- myadd_chunk(glb_chunks_df, "fit.data.training", major.inc=FALSE)
```

```
##                label step_major step_minor label_minor     bgn     end
## 14 fit.data.training          7          0           0 181.344 222.039
## 15 fit.data.training          7          1           1 222.039      NA
##    elapsed
## 14  40.695
## 15      NA
```


```r
#stop(here"); glb_to_sav()
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
                                    gsub(".", "\\.", glb_rsp_var_out, fixed = TRUE)),
                             "", mdlEnsembleComps)
    for (mdl_id in mdlEnsembleComps) {
        glbObsTrn <- glb_get_predictions(df=glbObsTrn, mdl_id=mdl_id, 
                                            rsp_var_out=glb_rsp_var_out,
                                            prob_threshold_def=prob_threshold)
        glbObsNew <- glb_get_predictions(df=glbObsNew, mdl_id=mdl_id, 
                                            rsp_var_out=glb_rsp_var_out,
                                            prob_threshold_def=prob_threshold)
    }    
}
glbObsTrn <- glb_get_predictions(df=glbObsTrn, mdl_id=glb_fin_mdl_id, 
                                     rsp_var_out=glb_rsp_var_out,
                                    prob_threshold_def=prob_threshold)
```

```
## Warning in glb_get_predictions(df = glbObsTrn, mdl_id = glb_fin_mdl_id, :
## Using default probability threshold: 0.1
```

```r
glb_featsimp_df <- myget_feats_importance(mdl=glb_fin_mdl,
                                          featsimp_df=glb_featsimp_df)
glb_featsimp_df[, paste0(glb_fin_mdl_id, ".importance")] <- glb_featsimp_df$importance
print(glb_featsimp_df)
```

```
##                                                    All.X.glmnet.importance
## NDSSName.my.fctrOpEd#Opinion#                                    100.00000
## NDSSName.my.fctrBusiness#Crosswords/Games#                        98.38173
## NDSSName.my.fctr#Opinion#ThePublicEditor                          87.34197
## NDSSName.my.fctrScience#Health#                                   84.05064
## NDSSName.my.fctrStyles#U.S.#                                      80.66804
## NDSSName.my.fctrBusiness#Technology#                              43.29623
## WordCount.log1p                                                   34.01597
## WordCount.root2                                                   32.29795
## .rnorm                                                            31.33944
## NDSSName.my.fctrCulture##                                         31.33944
## NDSSName.my.fctrCulture#Arts#                                     31.33944
## NDSSName.my.fctrMetro#N.Y./Region#                                31.33944
## WordCount.nexp                                                    31.33944
## NDSSName.my.fctrBusiness#BusinessDay#Dealbook                     31.33944
## NDSSName.my.fctrTravel#Travel#                                    31.31570
## NDSSName.my.fctrmyOther                                           27.85141
## NDSSName.my.fctrForeign#World#                                    30.00852
## NDSSName.my.fctr#Multimedia#                                      28.91940
## NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness                27.78546
## NDSSName.my.fctrStyles##Fashion                                   22.02991
## NDSSName.my.fctrTStyle##                                          17.36032
## NDSSName.my.fctrForeign#World#AsiaPacific                         19.29994
## NDSSName.my.fctr#U.S.#Education                                   18.25900
## NDSSName.my.fctr#Opinion#RoomForDebate                             0.00000
##                                                    importance
## NDSSName.my.fctrOpEd#Opinion#                       100.00000
## NDSSName.my.fctrBusiness#Crosswords/Games#           99.96082
## NDSSName.my.fctr#Opinion#ThePublicEditor             86.61487
## NDSSName.my.fctrScience#Health#                      82.33065
## NDSSName.my.fctrStyles#U.S.#                         77.12194
## NDSSName.my.fctrBusiness#Technology#                 38.62774
## WordCount.log1p                                      37.41260
## WordCount.root2                                      32.88366
## .rnorm                                               31.93272
## NDSSName.my.fctrCulture##                            31.93272
## NDSSName.my.fctrCulture#Arts#                        31.93272
## NDSSName.my.fctrMetro#N.Y./Region#                   31.93272
## WordCount.nexp                                       31.93272
## NDSSName.my.fctrBusiness#BusinessDay#Dealbook        31.12990
## NDSSName.my.fctrTravel#Travel#                       29.62778
## NDSSName.my.fctrmyOther                              26.20874
## NDSSName.my.fctrForeign#World#                       25.48953
## NDSSName.my.fctr#Multimedia#                         24.73162
## NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness   23.59037
## NDSSName.my.fctrStyles##Fashion                      21.10031
## NDSSName.my.fctrTStyle##                             16.67670
## NDSSName.my.fctrForeign#World#AsiaPacific            14.61070
## NDSSName.my.fctr#U.S.#Education                      14.34221
## NDSSName.my.fctr#Opinion#RoomForDebate                0.00000
##                                                    Final.glmnet.importance
## NDSSName.my.fctrOpEd#Opinion#                                    100.00000
## NDSSName.my.fctrBusiness#Crosswords/Games#                        99.96082
## NDSSName.my.fctr#Opinion#ThePublicEditor                          86.61487
## NDSSName.my.fctrScience#Health#                                   82.33065
## NDSSName.my.fctrStyles#U.S.#                                      77.12194
## NDSSName.my.fctrBusiness#Technology#                              38.62774
## WordCount.log1p                                                   37.41260
## WordCount.root2                                                   32.88366
## .rnorm                                                            31.93272
## NDSSName.my.fctrCulture##                                         31.93272
## NDSSName.my.fctrCulture#Arts#                                     31.93272
## NDSSName.my.fctrMetro#N.Y./Region#                                31.93272
## WordCount.nexp                                                    31.93272
## NDSSName.my.fctrBusiness#BusinessDay#Dealbook                     31.12990
## NDSSName.my.fctrTravel#Travel#                                    29.62778
## NDSSName.my.fctrmyOther                                           26.20874
## NDSSName.my.fctrForeign#World#                                    25.48953
## NDSSName.my.fctr#Multimedia#                                      24.73162
## NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness                23.59037
## NDSSName.my.fctrStyles##Fashion                                   21.10031
## NDSSName.my.fctrTStyle##                                          16.67670
## NDSSName.my.fctrForeign#World#AsiaPacific                         14.61070
## NDSSName.my.fctr#U.S.#Education                                   14.34221
## NDSSName.my.fctr#Opinion#RoomForDebate                             0.00000
```

```r
if (glb_is_classification && glb_is_binomial)
    glb_analytics_diag_plots(obs_df=glbObsTrn, mdl_id=glb_fin_mdl_id, 
            prob_threshold=glb_models_df[glb_models_df$id == glb_sel_mdl_id, 
                                         "opt.prob.threshold.OOB"]) else
    glb_analytics_diag_plots(obs_df=glbObsTrn, mdl_id=glb_fin_mdl_id)                  
```

![](NYTBlogs3_base_files/figure-html/fit.data.training_1-1.png) ![](NYTBlogs3_base_files/figure-html/fit.data.training_1-2.png) ![](NYTBlogs3_base_files/figure-html/fit.data.training_1-3.png) ![](NYTBlogs3_base_files/figure-html/fit.data.training_1-4.png) ![](NYTBlogs3_base_files/figure-html/fit.data.training_1-5.png) 

```
## [1] "Min/Max Boundaries: "
##      UniqueID Popular.fctr Popular.fctr.predict.All.X.glmnet.prob
## 96         96            N                                     NA
## 6370     6370            Y                              0.9608409
##      Popular.fctr.predict.All.X.glmnet
## 96                                <NA>
## 6370                                 Y
##      Popular.fctr.predict.All.X.glmnet.err
## 96                                      NA
## 6370                                 FALSE
##      Popular.fctr.predict.All.X.glmnet.err.abs
## 96                                          NA
## 6370                                0.03915913
##      Popular.fctr.predict.All.X.glmnet.accurate
## 96                                           NA
## 6370                                       TRUE
##      Popular.fctr.predict.Final.glmnet.prob
## 96                              0.006026625
## 6370                            0.915289043
##      Popular.fctr.predict.Final.glmnet
## 96                                   N
## 6370                                 Y
##      Popular.fctr.predict.Final.glmnet.err
## 96                                   FALSE
## 6370                                 FALSE
##      Popular.fctr.predict.Final.glmnet.err.abs
## 96                                 0.006026625
## 6370                               0.084710957
##      Popular.fctr.predict.Final.glmnet.accurate
## 96                                         TRUE
## 6370                                       TRUE
##      Popular.fctr.predict.Final.glmnet.error .label
## 96                                         0     96
## 6370                                       0   6370
## [1] "Inaccurate: "
##      UniqueID Popular.fctr Popular.fctr.predict.All.X.glmnet.prob
## 2182     2182            Y                            0.006633964
## 4020     4020            Y                                     NA
## 1696     1696            Y                            0.012803835
## 4775     4775            Y                                     NA
## 4721     4721            Y                            0.020278808
## 364       364            Y                            0.019247847
##      Popular.fctr.predict.All.X.glmnet
## 2182                                 N
## 4020                              <NA>
## 1696                                 N
## 4775                              <NA>
## 4721                                 N
## 364                                  N
##      Popular.fctr.predict.All.X.glmnet.err
## 2182                                  TRUE
## 4020                                    NA
## 1696                                  TRUE
## 4775                                    NA
## 4721                                  TRUE
## 364                                   TRUE
##      Popular.fctr.predict.All.X.glmnet.err.abs
## 2182                                 0.9933660
## 4020                                        NA
## 1696                                 0.9871962
## 4775                                        NA
## 4721                                 0.9797212
## 364                                  0.9807522
##      Popular.fctr.predict.All.X.glmnet.accurate
## 2182                                      FALSE
## 4020                                         NA
## 1696                                      FALSE
## 4775                                         NA
## 4721                                      FALSE
## 364                                       FALSE
##      Popular.fctr.predict.Final.glmnet.prob
## 2182                             0.01576970
## 4020                             0.02182292
## 1696                             0.02432199
## 4775                             0.03213481
## 4721                             0.03466036
## 364                              0.03629580
##      Popular.fctr.predict.Final.glmnet
## 2182                                 N
## 4020                                 N
## 1696                                 N
## 4775                                 N
## 4721                                 N
## 364                                  N
##      Popular.fctr.predict.Final.glmnet.err
## 2182                                  TRUE
## 4020                                  TRUE
## 1696                                  TRUE
## 4775                                  TRUE
## 4721                                  TRUE
## 364                                   TRUE
##      Popular.fctr.predict.Final.glmnet.err.abs
## 2182                                 0.9842303
## 4020                                 0.9781771
## 1696                                 0.9756780
## 4775                                 0.9678652
## 4721                                 0.9653396
## 364                                  0.9637042
##      Popular.fctr.predict.Final.glmnet.accurate
## 2182                                      FALSE
## 4020                                      FALSE
## 1696                                      FALSE
## 4775                                      FALSE
## 4721                                      FALSE
## 364                                       FALSE
##      Popular.fctr.predict.Final.glmnet.error
## 2182                             -0.08423030
## 4020                             -0.07817708
## 1696                             -0.07567801
## 4775                             -0.06786519
## 4721                             -0.06533964
## 364                              -0.06370420
##      UniqueID Popular.fctr Popular.fctr.predict.All.X.glmnet.prob
## 3029     3029            N                                     NA
## 323       323            N                             0.07847085
## 1181     1181            N                                     NA
## 2064     2064            N                             0.26251870
## 5287     5287            N                             0.60379585
## 4786     4786            N                             0.91965114
##      Popular.fctr.predict.All.X.glmnet
## 3029                              <NA>
## 323                                  N
## 1181                              <NA>
## 2064                                 Y
## 5287                                 Y
## 4786                                 Y
##      Popular.fctr.predict.All.X.glmnet.err
## 3029                                    NA
## 323                                  FALSE
## 1181                                    NA
## 2064                                  TRUE
## 5287                                  TRUE
## 4786                                  TRUE
##      Popular.fctr.predict.All.X.glmnet.err.abs
## 3029                                        NA
## 323                                 0.07847085
## 1181                                        NA
## 2064                                0.26251870
## 5287                                0.60379585
## 4786                                0.91965114
##      Popular.fctr.predict.All.X.glmnet.accurate
## 3029                                         NA
## 323                                        TRUE
## 1181                                         NA
## 2064                                      FALSE
## 5287                                      FALSE
## 4786                                      FALSE
##      Popular.fctr.predict.Final.glmnet.prob
## 3029                              0.1051303
## 323                               0.1168606
## 1181                              0.1702228
## 2064                              0.2950796
## 5287                              0.4139998
## 4786                              0.8143734
##      Popular.fctr.predict.Final.glmnet
## 3029                                 Y
## 323                                  Y
## 1181                                 Y
## 2064                                 Y
## 5287                                 Y
## 4786                                 Y
##      Popular.fctr.predict.Final.glmnet.err
## 3029                                  TRUE
## 323                                   TRUE
## 1181                                  TRUE
## 2064                                  TRUE
## 5287                                  TRUE
## 4786                                  TRUE
##      Popular.fctr.predict.Final.glmnet.err.abs
## 3029                                 0.1051303
## 323                                  0.1168606
## 1181                                 0.1702228
## 2064                                 0.2950796
## 5287                                 0.4139998
## 4786                                 0.8143734
##      Popular.fctr.predict.Final.glmnet.accurate
## 3029                                      FALSE
## 323                                       FALSE
## 1181                                      FALSE
## 2064                                      FALSE
## 5287                                      FALSE
## 4786                                      FALSE
##      Popular.fctr.predict.Final.glmnet.error
## 3029                             0.005130305
## 323                              0.016860557
## 1181                             0.070222845
## 2064                             0.195079625
## 5287                             0.313999830
## 4786                             0.714373384
##      UniqueID Popular.fctr Popular.fctr.predict.All.X.glmnet.prob
## 221       221            N                                     NA
## 3590     3590            N                                     NA
## 472       472            N                                     NA
## 2995     2995            N                                     NA
## 6276     6276            N                                     NA
## 3258     3258            N                              0.9735424
##      Popular.fctr.predict.All.X.glmnet
## 221                               <NA>
## 3590                              <NA>
## 472                               <NA>
## 2995                              <NA>
## 6276                              <NA>
## 3258                                 Y
##      Popular.fctr.predict.All.X.glmnet.err
## 221                                     NA
## 3590                                    NA
## 472                                     NA
## 2995                                    NA
## 6276                                    NA
## 3258                                  TRUE
##      Popular.fctr.predict.All.X.glmnet.err.abs
## 221                                         NA
## 3590                                        NA
## 472                                         NA
## 2995                                        NA
## 6276                                        NA
## 3258                                 0.9735424
##      Popular.fctr.predict.All.X.glmnet.accurate
## 221                                          NA
## 3590                                         NA
## 472                                          NA
## 2995                                         NA
## 6276                                         NA
## 3258                                      FALSE
##      Popular.fctr.predict.Final.glmnet.prob
## 221                               0.9169618
## 3590                              0.9171703
## 472                               0.9195159
## 2995                              0.9203097
## 6276                              0.9215804
## 3258                              0.9245769
##      Popular.fctr.predict.Final.glmnet
## 221                                  Y
## 3590                                 Y
## 472                                  Y
## 2995                                 Y
## 6276                                 Y
## 3258                                 Y
##      Popular.fctr.predict.Final.glmnet.err
## 221                                   TRUE
## 3590                                  TRUE
## 472                                   TRUE
## 2995                                  TRUE
## 6276                                  TRUE
## 3258                                  TRUE
##      Popular.fctr.predict.Final.glmnet.err.abs
## 221                                  0.9169618
## 3590                                 0.9171703
## 472                                  0.9195159
## 2995                                 0.9203097
## 6276                                 0.9215804
## 3258                                 0.9245769
##      Popular.fctr.predict.Final.glmnet.accurate
## 221                                       FALSE
## 3590                                      FALSE
## 472                                       FALSE
## 2995                                      FALSE
## 6276                                      FALSE
## 3258                                      FALSE
##      Popular.fctr.predict.Final.glmnet.error
## 221                                0.8169618
## 3590                               0.8171703
## 472                                0.8195159
## 2995                               0.8203097
## 6276                               0.8215804
## 3258                               0.8245769
```

![](NYTBlogs3_base_files/figure-html/fit.data.training_1-6.png) 

```r
dsp_feats_vctr <- c(NULL)
for(var in grep(".importance", names(glb_feats_df), fixed=TRUE, value=TRUE))
    dsp_feats_vctr <- union(dsp_feats_vctr, 
                            glb_feats_df[!is.na(glb_feats_df[, var]), "id"])

# print(glbObsTrn[glbObsTrn$UniqueID %in% FN_OOB_ids, 
#                     grep(glb_rsp_var, names(glbObsTrn), value=TRUE)])

print(setdiff(names(glbObsTrn), names(glbObsAll)))
```

```
## [1] "Popular.fctr.predict.Final.glmnet.prob"    
## [2] "Popular.fctr.predict.Final.glmnet"         
## [3] "Popular.fctr.predict.Final.glmnet.err"     
## [4] "Popular.fctr.predict.Final.glmnet.err.abs" 
## [5] "Popular.fctr.predict.Final.glmnet.accurate"
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

![](NYTBlogs3_base_files/figure-html/fit.data.training_1-7.png) 

```r
glb_chunks_df <- myadd_chunk(glb_chunks_df, "predict.data.new", major.inc=TRUE)
```

```
##                label step_major step_minor label_minor     bgn     end
## 15 fit.data.training          7          1           1 222.039 233.199
## 16  predict.data.new          8          0           0 233.199      NA
##    elapsed
## 15   11.16
## 16      NA
```

## Step `8.0: predict data new`

```r
# Compute final model predictions

#glb_to_sav(); all.equal(sav_allobs_df, glbObsAll); all.equal(sav_trnobs_df, glbObsTrn); all.equal(sav_newobs_df, glbObsNew)  
if (glb_is_classification && glb_is_binomial)
    prob_threshold_def <- 
        glb_models_df[glb_models_df$id == glb_sel_mdl_id, "opt.prob.threshold.OOB"] else
    prob_threshold_def <- NULL
for (obsSet in c("trn", "new")) {
    obs_df <- switch(obsSet, all = glbObsAll, trn = glbObsTrn, new = glbObsNew)
    obs_df <- glb_get_predictions(obs_df, mdl_id = glb_fin_mdl_id, 
                    rsp_var_out = glb_rsp_var_out, prob_threshold_def = prob_threshold_def)
    if (obsSet == "all") glbObsAll <- obs_df else
    if (obsSet == "trn") glbObsTrn <- obs_df else
    if (obsSet == "new") glbObsNew <- obs_df
}
```

```
## Warning in glb_get_predictions(obs_df, mdl_id = glb_fin_mdl_id, rsp_var_out
## = glb_rsp_var_out, : Using default probability threshold: 0.1
```

```
## Warning in glb_get_predictions(obs_df, mdl_id = glb_fin_mdl_id, rsp_var_out
## = glb_rsp_var_out, : Using default probability threshold: 0.1
```

```r
rm(obs_df)
glbObsAll <- orderBy(reformulate(glb_id_var), myrbind_df(glbObsTrn, glbObsNew))

glb_analytics_diag_plots(obs_df = glbObsNew, mdl_id = glb_fin_mdl_id, 
                         prob_threshold = prob_threshold_def)
```

```
## Warning: Removed 1870 rows containing missing values (geom_point).
```

```
## Warning: Removed 1870 rows containing missing values (geom_point).
```

![](NYTBlogs3_base_files/figure-html/predict.data.new-1.png) 

```
## Warning: Removed 1870 rows containing missing values (geom_point).
```

```
## Warning: Removed 1870 rows containing missing values (geom_point).
```

![](NYTBlogs3_base_files/figure-html/predict.data.new-2.png) 

```
## Warning: Removed 1870 rows containing missing values (geom_point).
```

```
## Warning: Removed 1870 rows containing missing values (geom_point).
```

![](NYTBlogs3_base_files/figure-html/predict.data.new-3.png) 

```
## Warning: Removed 1870 rows containing missing values (geom_point).
```

```
## Warning: Removed 1870 rows containing missing values (geom_point).
```

![](NYTBlogs3_base_files/figure-html/predict.data.new-4.png) 

```
## Warning: Removed 1870 rows containing missing values (geom_point).
```

```
## Warning: Removed 1870 rows containing missing values (geom_point).
```

![](NYTBlogs3_base_files/figure-html/predict.data.new-5.png) 

```
## NULL
```

```r
if (is.null(glb_out_obs)) obs_df <- glbObsNew else
    obs_df <- switch(glb_out_obs, 
                     all = glbObsAll, trn = glbObsTrn, new = glbObsNew)

require(stringr)
```

```
## Loading required package: stringr
```

```r
if (glb_is_classification && glb_is_binomial) {
    obsout_df <- obs_df[, glb_id_var, FALSE]
    for (clmn in names(glb_out_vars_lst))
        if (!grepl("^%<d-%", glb_out_vars_lst[[clmn]]))
            obsout_df[, clmn] <- obs_df[, glb_out_vars_lst[[clmn]]] else {
            feat <- str_trim(unlist(strsplit(glb_out_vars_lst[[clmn]], "%<d-%"))[2])
            obsout_df[, clmn] <- obs_df[, eval(parse(text = feat))]
        }                                        
    
#     glb_force_prediction_lst <- list()
#     glb_force_prediction_lst[["0"]] <- c(11885, 11907, 11932, 11943, 
#                                          12050, 12115, 12171, 
#                                          12253, 12285, 12367, 12388, 12399,
#                                          12585)
#     for (obs_id in glb_force_prediction_lst[["0"]]) {
#         if (sum(glbObsAll[, glb_id_var] == obs_id) == 0)
#             next
#         if (is.na(glbObsAll[glbObsAll[, glb_id_var] == obs_id, ".grpid"]))
#             stop(".grpid is NA")
# #         submit_df[submit_df[, glb_id_var] == obs_id, "Probability1"] <-
# #             max(0, submit_df[submit_df[, glb_id_var] == obs_id, "Probability1"] - 0.5)
#     }    
    
#     glb_force_prediction_lst[["1"]] <- c(11871, 11875, 11886, 
#                     11913, 11931, 11937, 11967, 11982, 11990, 11991, 11994, 11999,
#                                 12000, 12002, 12004, 12018, 12021, 12065, 12072,
#                                          12111, 12114, 12126, 12134, 12152, 12172,
#                                          12213, 12214, 12233, 12265, 12278, 12299, 
#                                          12446, 12491, 
#                                          12505, 12576, 12608, 12630)
#     for (obs_id in glb_force_prediction_lst[["1"]]) {
#         if (sum(glbObsAll[, glb_id_var] == obs_id) == 0)
#             next
#         if (is.na(glbObsAll[glbObsAll[, glb_id_var] == obs_id, ".grpid"]))
#             stop(".grpid is NA")
# #         submit_df[submit_df[, glb_id_var] == obs_id, "Probability1"] <-
# #             min(0.9999, submit_df[submit_df[, glb_id_var] == obs_id, "Probability1"] + 0.5)
#     }    
    
#     rsp_var_out <- paste0(glb_rsp_var_out, glb_fin_mdl_id)
#     for (obs_id in glbObsNew[!is.na(glbObsNew[, rsp_var_out]) & 
#                                  (glbObsNew[, rsp_var_out] == "Y") & 
#                                  (glbObsNew[ , "startprice"] > 675), "UniqueID"]) {
# #         submit_df[submit_df[, glb_id_var] == obs_id, "Probability1"] <-
# #             max(0, submit_df[submit_df[, glb_id_var] == obs_id, "Probability1"] - 0.5)
#     }    
} else {
#     submit_df <- glbObsNew[, c(glb_id_var, 
#                                    paste0(glb_rsp_var_out, glb_fin_mdl_id))]
    obsout_df <- obs_df[, glb_id_var, FALSE]
    for (clmn in names(glb_out_vars_lst))
        if (!grepl("^%<d-%", glb_out_vars_lst[[clmn]]))
            obsout_df[, clmn] <- obs_df[, glb_out_vars_lst[[clmn]]] else {
            feat <- str_trim(unlist(strsplit(glb_out_vars_lst[[clmn]], "%<d-%"))[2])
            obsout_df[, clmn] <- obs_df[, eval(parse(text=feat))]
        }                                        
}    

if (glb_is_classification) {
    rsp_var_out <- paste0(glb_rsp_var_out, glb_fin_mdl_id)
    if (".grpid" %in% names(glbObsNew)) {
        # Dups were found in glbObsAll
        tmp_newobs_df <- subset(glbObsNew[, c(glb_id_var, ".grpid", rsp_var_out)],
                                !is.na(.grpid))
        tmp_newobs_df <- 
            merge(tmp_newobs_df, dupgrps_df, by = ".grpid", all.x = TRUE)
        tmp_newobs_df <- 
            merge(tmp_newobs_df, obsout_df, by = glb_id_var, all.x = TRUE)
        tmp_newobs_df$.err <- 
            ((tmp_newobs_df$Probability1 > 0.5) & (tmp_newobs_df$sold.0 > 0) |
             (tmp_newobs_df$Probability1 < 0.5) & (tmp_newobs_df$sold.1 > 0))
        tmp_newobs_df <- orderBy(~UniqueID, subset(tmp_newobs_df, .err == TRUE))
        print(sprintf("Prediction errors in duplicates: %d", nrow(tmp_newobs_df)))
        print(tmp_newobs_df)
    }
    
    tmp_newobs_df <- cbind(glbObsNew, obsout_df[, "Probability1", FALSE])

    # Check predictions that are outside of data ranges
#stop(here")    
    require(stringr)
    tmp_feats_df <- subset(glb_feats_df, 
                           !nzv & 
                            (exclude.as.feat != 1) & 
                            !grepl(".fctr", id, fixed=TRUE))[, "id", FALSE]
    ranges_all_df <- glbObsAll[, tmp_feats_df$id] %>% 
                        dplyr::summarise_each(funs(min(., na.rm=TRUE), 
                                                   max(., na.rm=TRUE))) %>%
                        tidyr::gather() %>%
                        dplyr::mutate(id=str_sub(key, 1, -5), 
                                      stat=str_sub(key, -3)) %>% 
                        dplyr::select(-key) %>%
                        tidyr::spread(stat, value)
    
#     sav_ranges_trn_df <- ranges_trn_df; all.equal(sav_ranges_trn_df, ranges_trn_df)
#     sav_ranges_new_df <- ranges_new_df; all.equal(sav_ranges_new_df, ranges_new_df)    
    get_ranges_df <- function(obs_df, feats, class_var) {
        require(tidyr)
        ranges_df <- obs_df[, c(class_var, feats)] %>% 
            dplyr::group_by_(class_var) %>%
            dplyr::summarise_each(funs(min(., na.rm=TRUE), 
                                       max(., na.rm=TRUE))) %>%
            tidyr::gather(key, value, -1) %>%
            mutate(id=str_sub(key, 1, -5), 
                   stat.vname=paste0(str_sub(key, -3), ".", class_var)) %>%
            unite_("stat.class", c("stat.vname", class_var), sep=".") %>% 
            dplyr::select(-key) %>%
            spread(stat.class, value)
        return(ranges_df)
    }
    rsp_var_out_OOB <- paste0(glb_rsp_var_out, glb_sel_mdl_id)
    rsp_var_out_new <- paste0(glb_rsp_var_out, glb_fin_mdl_id)    
    ranges_trn_df <- get_ranges_df(obs_df=glbObsTrn, feats=tmp_feats_df$id, 
                                   class_var=glb_rsp_var)
    ranges_fit_df <- get_ranges_df(obs_df=glbObsFit, feats=tmp_feats_df$id, 
                                   class_var=glb_rsp_var)
    ranges_OOB_df <- get_ranges_df(obs_df=glbObsOOB, feats=tmp_feats_df$id, 
                                   class_var=rsp_var_out_OOB)
    ranges_new_df <- get_ranges_df(obs_df=glbObsNew, feats=tmp_feats_df$id, 
                                   class_var=rsp_var_out_new)

    for (obsset in c("OOB", "new")) {
        if (obsset == "OOB") { 
            ranges_ref_df <- ranges_fit_df; obs_df <- glbObsOOB; 
            rsp_var_out_obs <- rsp_var_out_OOB; sprintf_pfx <- "OOBobs";
        } else { 
            ranges_ref_df <- ranges_trn_df; obs_df <- glbObsNew; 
            rsp_var_out_obs <- rsp_var_out_new; sprintf_pfx <- "newobs"; 
        }
        plt_feats_df <- glb_feats_df %>% 
                            merge(ranges_all_df, all=TRUE) %>%
                            merge(ranges_ref_df, all=TRUE) %>%
                            merge(ranges_OOB_df, all=TRUE) %>%        
                            merge(ranges_new_df, all=TRUE) %>%
                            subset(!is.na(min) & (id != ".rnorm"))
        row.names(plt_feats_df) <- plt_feats_df$id
        range_outlier_ids <- c(NULL)
        for (clss in unique(obs_df[, rsp_var_out_obs])) {
            for (stat in c("min", "max")) {
                if (stat == "min") {
                    dsp_feats <- plt_feats_df[
                            which(plt_feats_df[, paste("min", rsp_var_out_obs, clss, sep=".")] < 
                                  plt_feats_df[, paste("min", glb_rsp_var, clss, sep=".")]), "id"]
                } else {
                    dsp_feats <- plt_feats_df[
                            which(plt_feats_df[, paste("max", rsp_var_out_obs, clss, sep=".")] > 
                                  plt_feats_df[, paste("max", glb_rsp_var, clss, sep=".")]), "id"]
                }
                if (length(dsp_feats) > 0) {
                    ths_ids <- c(NULL)
                    for (feat in dsp_feats) {
                        if (stat == "min") {
                            ths_ids <- union(ths_ids, 
                                             obs_df[(obs_df[, rsp_var_out_obs] == clss) &
                                                           (obs_df[, feat] < 
                plt_feats_df[plt_feats_df$id == feat, paste("min", glb_rsp_var, clss, sep=".")]), 
                                                            glb_id_var])
                        } else {
                        ths_ids <- union(ths_ids, 
                                             obs_df[(obs_df[, rsp_var_out_obs] == clss) &
                                                           (obs_df[, feat] > 
                plt_feats_df[plt_feats_df$id == feat, paste("max", glb_rsp_var, clss, sep=".")]), 
                                                            glb_id_var])
                        }
                    }
                    tmp_obs_df <- obs_df[obs_df[, glb_id_var] %in% ths_ids, 
                                                   c(glb_id_var, rsp_var_out_obs, dsp_feats)]
                    if (stat == "min") {
                        print(sprintf("%s %s %s: min < min of Train range: %d", 
                                      sprintf_pfx, rsp_var_out_obs, clss, nrow(tmp_obs_df)))
                    } else {
                        print(sprintf("%s %s %s: max > max of Train range: %d", 
                                      sprintf_pfx, rsp_var_out_obs, clss, nrow(tmp_obs_df)))
                    }
                    myprint_df(tmp_obs_df)
                    print(subset(plt_feats_df, id %in% dsp_feats))
                    
                    range_outlier_ids <- union(range_outlier_ids, ths_ids)
                }
            }
        }
        print(sprintf("%s total range outliers: %d", sprintf_pfx,
                      length(range_outlier_ids)))
    }
}
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
## [1] "OOBobs Popular.fctr.predict.All.X.glmnet Y: min < min of Train range: 1"
##      UniqueID Popular.fctr.predict.All.X.glmnet WordCount.log1p
## 6435     6435                                 Y               0
##      WordCount.root2
## 6435               0
##                              id     cor.y exclude.as.feat cor.y.abs
## WordCount.log1p WordCount.log1p 0.2543196           FALSE 0.2543196
## WordCount.root2 WordCount.root2 0.2921207           FALSE 0.2921207
##                      cor.high.X freqRatio percentUnique zeroVar   nzv
## WordCount.log1p WordCount.root2  2.315789      24.15799   FALSE FALSE
## WordCount.root2            <NA>  2.315789      24.15799   FALSE FALSE
##                 is.cor.y.abs.low interaction.feat shapiro.test.p.value
## WordCount.log1p            FALSE               NA         1.576866e-49
## WordCount.root2            FALSE               NA         4.556481e-30
##                 rsp_var_raw id_var rsp_var       max min
## WordCount.log1p       FALSE     NA      NA   9.29771   0
## WordCount.root2       FALSE     NA      NA 104.46052   0
##                 max.Popular.fctr.N max.Popular.fctr.Y min.Popular.fctr.N
## WordCount.log1p           8.819665            9.29771                  0
## WordCount.root2          82.249620          104.46052                  0
##                 min.Popular.fctr.Y max.Popular.fctr.predict.All.X.glmnet.N
## WordCount.log1p            1.94591                                8.117014
## WordCount.root2            2.44949                               57.879185
##                 max.Popular.fctr.predict.All.X.glmnet.Y
## WordCount.log1p                                9.140883
## WordCount.root2                               96.581572
##                 min.Popular.fctr.predict.All.X.glmnet.N
## WordCount.log1p                                       0
## WordCount.root2                                       0
##                 min.Popular.fctr.predict.All.X.glmnet.Y
## WordCount.log1p                                       0
## WordCount.root2                                       0
##                 max.Popular.fctr.predict.Final.glmnet.N
## WordCount.log1p                                7.797702
## WordCount.root2                               49.335586
##                 max.Popular.fctr.predict.Final.glmnet.Y
## WordCount.log1p                                8.692322
## WordCount.root2                               77.175126
##                 min.Popular.fctr.predict.Final.glmnet.N
## WordCount.log1p                                       0
## WordCount.root2                                       0
##                 min.Popular.fctr.predict.Final.glmnet.Y
## WordCount.log1p                                1.609438
## WordCount.root2                                2.000000
## [1] "OOBobs Popular.fctr.predict.All.X.glmnet Y: max > max of Train range: 1"
##      UniqueID Popular.fctr.predict.All.X.glmnet WordCount.nexp
## 6435     6435                                 Y              1
##                            id      cor.y exclude.as.feat cor.y.abs
## WordCount.nexp WordCount.nexp -0.0532084           FALSE 0.0532084
##                cor.high.X freqRatio percentUnique zeroVar   nzv
## WordCount.nexp       <NA>  17.76136      11.32884   FALSE FALSE
##                is.cor.y.abs.low interaction.feat shapiro.test.p.value
## WordCount.nexp            FALSE               NA         9.108805e-94
##                rsp_var_raw id_var rsp_var max min max.Popular.fctr.N
## WordCount.nexp       FALSE     NA      NA   1   0                  1
##                max.Popular.fctr.Y min.Popular.fctr.N min.Popular.fctr.Y
## WordCount.nexp        0.002478752                  0                  0
##                max.Popular.fctr.predict.All.X.glmnet.N
## WordCount.nexp                                       1
##                max.Popular.fctr.predict.All.X.glmnet.Y
## WordCount.nexp                                       1
##                min.Popular.fctr.predict.All.X.glmnet.N
## WordCount.nexp                                       0
##                min.Popular.fctr.predict.All.X.glmnet.Y
## WordCount.nexp                                       0
##                max.Popular.fctr.predict.Final.glmnet.N
## WordCount.nexp                                       1
##                max.Popular.fctr.predict.Final.glmnet.Y
## WordCount.nexp                              0.01831564
##                min.Popular.fctr.predict.Final.glmnet.N
## WordCount.nexp                                       0
##                min.Popular.fctr.predict.Final.glmnet.Y
## WordCount.nexp                                       0
## [1] "OOBobs total range outliers: 1"
## [1] "newobs Popular.fctr.predict.Final.glmnet Y: min < min of Train range: 1"
##      UniqueID Popular.fctr.predict.Final.glmnet WordCount.log1p
## 8217     8217                                 Y        1.609438
##      WordCount.root2
## 8217               2
##                              id     cor.y exclude.as.feat cor.y.abs
## WordCount.log1p WordCount.log1p 0.2543196           FALSE 0.2543196
## WordCount.root2 WordCount.root2 0.2921207           FALSE 0.2921207
##                      cor.high.X freqRatio percentUnique zeroVar   nzv
## WordCount.log1p WordCount.root2  2.315789      24.15799   FALSE FALSE
## WordCount.root2            <NA>  2.315789      24.15799   FALSE FALSE
##                 is.cor.y.abs.low interaction.feat shapiro.test.p.value
## WordCount.log1p            FALSE               NA         1.576866e-49
## WordCount.root2            FALSE               NA         4.556481e-30
##                 rsp_var_raw id_var rsp_var       max min
## WordCount.log1p       FALSE     NA      NA   9.29771   0
## WordCount.root2       FALSE     NA      NA 104.46052   0
##                 max.Popular.fctr.N max.Popular.fctr.Y min.Popular.fctr.N
## WordCount.log1p           8.819665            9.29771                  0
## WordCount.root2          82.249620          104.46052                  0
##                 min.Popular.fctr.Y max.Popular.fctr.predict.All.X.glmnet.N
## WordCount.log1p            1.94591                                8.117014
## WordCount.root2            2.44949                               57.879185
##                 max.Popular.fctr.predict.All.X.glmnet.Y
## WordCount.log1p                                9.140883
## WordCount.root2                               96.581572
##                 min.Popular.fctr.predict.All.X.glmnet.N
## WordCount.log1p                                       0
## WordCount.root2                                       0
##                 min.Popular.fctr.predict.All.X.glmnet.Y
## WordCount.log1p                                       0
## WordCount.root2                                       0
##                 max.Popular.fctr.predict.Final.glmnet.N
## WordCount.log1p                                7.797702
## WordCount.root2                               49.335586
##                 max.Popular.fctr.predict.Final.glmnet.Y
## WordCount.log1p                                8.692322
## WordCount.root2                               77.175126
##                 min.Popular.fctr.predict.Final.glmnet.N
## WordCount.log1p                                       0
## WordCount.root2                                       0
##                 min.Popular.fctr.predict.Final.glmnet.Y
## WordCount.log1p                                1.609438
## WordCount.root2                                2.000000
## [1] "newobs Popular.fctr.predict.Final.glmnet Y: max > max of Train range: 1"
##      UniqueID Popular.fctr.predict.Final.glmnet WordCount.nexp
## 8217     8217                                 Y     0.01831564
##                            id      cor.y exclude.as.feat cor.y.abs
## WordCount.nexp WordCount.nexp -0.0532084           FALSE 0.0532084
##                cor.high.X freqRatio percentUnique zeroVar   nzv
## WordCount.nexp       <NA>  17.76136      11.32884   FALSE FALSE
##                is.cor.y.abs.low interaction.feat shapiro.test.p.value
## WordCount.nexp            FALSE               NA         9.108805e-94
##                rsp_var_raw id_var rsp_var max min max.Popular.fctr.N
## WordCount.nexp       FALSE     NA      NA   1   0                  1
##                max.Popular.fctr.Y min.Popular.fctr.N min.Popular.fctr.Y
## WordCount.nexp        0.002478752                  0                  0
##                max.Popular.fctr.predict.All.X.glmnet.N
## WordCount.nexp                                       1
##                max.Popular.fctr.predict.All.X.glmnet.Y
## WordCount.nexp                                       1
##                min.Popular.fctr.predict.All.X.glmnet.N
## WordCount.nexp                                       0
##                min.Popular.fctr.predict.All.X.glmnet.Y
## WordCount.nexp                                       0
##                max.Popular.fctr.predict.Final.glmnet.N
## WordCount.nexp                                       1
##                max.Popular.fctr.predict.Final.glmnet.Y
## WordCount.nexp                              0.01831564
##                min.Popular.fctr.predict.Final.glmnet.N
## WordCount.nexp                                       0
##                min.Popular.fctr.predict.Final.glmnet.Y
## WordCount.nexp                                       0
## [1] "newobs total range outliers: 1"
```

```r
#stop(here"); glb_to_sav(); sav_obsout_df <- obsout_df; all.equal(sav_obsout_df, obsout_df); obsout_df <- sav_obsout_df

# This does not work for classification since AUC distribution might be different for different models
#   -> Run glm on .prob from this glb_fin_mdl_id & .prob from stacked file & stack condition as a feature
if (!is.null(glbOutStackFnames)) {
    for (fname in glbOutStackFnames) {
        print(sprintf("Stacking file %s to prediction output...", fname))
        #obsout_df <- dplyr::arrange_(rbind(obsout_df, read.csv(fname)), "UniqueID")
        obsout_df <- dplyr::arrange_(rbind(obsout_df, 
                #read.csv(fname) %>% filter(!(UniqueID %in% obsout_df$UniqueID))),
            #read.csv(fname) %>% filter(!(UniqueID %in% obsout_df[, glb_id_var]))),
            read.csv(fname) %>% 
                dplyr::filter_(interp(~!(var %in% obsout_df$var), 
                                      var = as.name(glb_id_var)))),
                                    glb_id_var)
        
        if (nrow(obsout_df) != length(unique(obsout_df[, glb_id_var])))
            stop("Potential dups in stacked prediction output")
    }
}

out_fname <- paste0(glb_out_pfx, "out.csv")
write.csv(obsout_df, out_fname, quote = FALSE, row.names = FALSE)
#cat(" ", "\n", file=submit_fn, append=TRUE)

# print(orderBy(~ -max.auc.OOB, glb_models_df[, c("model_id", 
#             "max.auc.OOB", "max.Accuracy.OOB")]))
for (txt_var in glbFeatsText) {
    # Print post-stem-words but need post-stop-words for debugging ?
    print(sprintf("    All post-stem-words TfIDf terms for %s:", txt_var))
    myprint_df(glb_post_stem_words_terms_df_lst[[txt_var]])
    TfIdf_mtrx <- glb_post_stem_words_TfIdf_mtrx_lst[[txt_var]]
    print(glbObsAll[
        which(TfIdf_mtrx[, tail(glb_post_stem_words_terms_df_lst[[txt_var]], 1)$pos] > 0), 
                        c(glb_id_var, glbFeatsText)])
    print(nrow(subset(glb_post_stem_words_terms_df_lst[[txt_var]], freq == 1)))
    #print(glbObsAll[which(TfIdf_mtrx[, 207] > 0), c(glb_id_var, glbFeatsText)])
    #unlist(strsplit(glbObsAll[2157, "description"], ""))
    #glbObsAll[2442, c(glb_id_var, glbFeatsText)]
    #TfIdf_mtrx[2442, TfIdf_mtrx[2442, ] > 0]  

    print(sprintf("    Top_n post_stem_words TfIDf terms for %s:", txt_var))
    tmp_df <- glb_post_stem_words_terms_df_lst[[txt_var]]
    top_n_vctr <- tmp_df$term[1:glb_txt_top_n[[txt_var]]]
    tmp_freq1_df <- subset(tmp_df, freq == 1)
    tmp_freq1_df$top_n <- grepl(paste0(top_n_vctr, collapse="|"), tmp_freq1_df$term)
    print(subset(tmp_freq1_df, top_n == TRUE))
}

if (glb_is_classification && glb_is_binomial)
    print(glb_models_df[glb_models_df$model_id == glb_sel_mdl_id, 
                        "opt.prob.threshold.OOB"])
```

```
## numeric(0)
```

```r
print(sprintf("glb_sel_mdl_id: %s", glb_sel_mdl_id))
```

```
## [1] "glb_sel_mdl_id: All.X.glmnet"
```

```r
print(sprintf("glb_fin_mdl_id: %s", glb_fin_mdl_id))
```

```
## [1] "glb_fin_mdl_id: Final.glmnet"
```

```r
print(dsp_models_df)
```

```
##                                                        id max.Accuracy.OOB
## Max.cor.Y.rpart                           Max.cor.Y.rpart        0.8200231
## All.X.glmnet                                 All.X.glmnet        0.7690972
## Max.cor.Y.rcv.1X1.cp.0.rpart Max.cor.Y.rcv.1X1.cp.0.rpart        0.7673611
## Max.cor.Y.rcv.1X1.glmnet         Max.cor.Y.rcv.1X1.glmnet        0.7604167
## Max.cor.Y.rcv.5X3.glmnet         Max.cor.Y.rcv.5X3.glmnet        0.7604167
## Max.cor.Y.rcv.5X1.glmnet         Max.cor.Y.rcv.5X1.glmnet        0.7604167
## Max.cor.Y.rcv.5X5.glmnet         Max.cor.Y.rcv.5X5.glmnet        0.7604167
## Max.cor.Y.rcv.3X1.glmnet         Max.cor.Y.rcv.3X1.glmnet        0.7575231
## Max.cor.Y.rcv.3X3.glmnet         Max.cor.Y.rcv.3X3.glmnet        0.7575231
## Interact.High.cor.Y.glmnet     Interact.High.cor.Y.glmnet        0.7575231
## Low.cor.X.glmnet                         Low.cor.X.glmnet        0.7575231
## Max.cor.Y.rcv.3X5.glmnet         Max.cor.Y.rcv.3X5.glmnet        0.7575231
## MFO.myMFO_classfr                       MFO.myMFO_classfr        0.1331019
## Random.myrandom_classfr           Random.myrandom_classfr        0.1331019
##                              max.AUCROCR.OOB max.AUCpROC.OOB
## Max.cor.Y.rpart                    0.5892132       0.5870523
## All.X.glmnet                       0.8061009       0.5965780
## Max.cor.Y.rcv.1X1.cp.0.rpart       0.7773858       0.6174697
## Max.cor.Y.rcv.1X1.glmnet           0.8116126       0.5962443
## Max.cor.Y.rcv.5X3.glmnet           0.8114863       0.5962443
## Max.cor.Y.rcv.5X1.glmnet           0.8114863       0.5962443
## Max.cor.Y.rcv.5X5.glmnet           0.8114863       0.5962443
## Max.cor.Y.rcv.3X1.glmnet           0.8067975       0.5962443
## Max.cor.Y.rcv.3X3.glmnet           0.8067975       0.5962443
## Interact.High.cor.Y.glmnet         0.8067975       0.5962443
## Low.cor.X.glmnet                   0.8067975       0.5962443
## Max.cor.Y.rcv.3X5.glmnet           0.8067975       0.5962443
## MFO.myMFO_classfr                  0.5000000       0.5000000
## Random.myrandom_classfr            0.4857956       0.5125675
##                              max.Accuracy.fit opt.prob.threshold.fit
## Max.cor.Y.rpart                     0.9296422                    0.6
## All.X.glmnet                        0.9326952                    0.4
## Max.cor.Y.rcv.1X1.cp.0.rpart        0.9381765                    0.4
## Max.cor.Y.rcv.1X1.glmnet            0.9329725                    0.5
## Max.cor.Y.rcv.5X3.glmnet            0.9333905                    0.5
## Max.cor.Y.rcv.5X1.glmnet            0.9331818                    0.5
## Max.cor.Y.rcv.5X5.glmnet            0.9331816                    0.5
## Max.cor.Y.rcv.3X1.glmnet            0.9335973                    0.4
## Max.cor.Y.rcv.3X3.glmnet            0.9333193                    0.4
## Interact.High.cor.Y.glmnet          0.9333193                    0.4
## Low.cor.X.glmnet                    0.9333193                    0.4
## Max.cor.Y.rcv.3X5.glmnet            0.9332218                    0.4
## MFO.myMFO_classfr                   0.1796420                    0.1
## Random.myrandom_classfr             0.1796420                    0.1
##                              opt.prob.threshold.OOB
## Max.cor.Y.rpart                                 0.6
## All.X.glmnet                                    0.1
## Max.cor.Y.rcv.1X1.cp.0.rpart                    0.1
## Max.cor.Y.rcv.1X1.glmnet                        0.1
## Max.cor.Y.rcv.5X3.glmnet                        0.1
## Max.cor.Y.rcv.5X1.glmnet                        0.1
## Max.cor.Y.rcv.5X5.glmnet                        0.1
## Max.cor.Y.rcv.3X1.glmnet                        0.1
## Max.cor.Y.rcv.3X3.glmnet                        0.1
## Interact.High.cor.Y.glmnet                      0.1
## Low.cor.X.glmnet                                0.1
## Max.cor.Y.rcv.3X5.glmnet                        0.1
## MFO.myMFO_classfr                               0.1
## Random.myrandom_classfr                         0.1
```

```r
if (glb_is_regression) {
    print(sprintf("%s OOB RMSE: %0.4f", glb_sel_mdl_id,
            glb_models_df[glb_models_df$model_id == glb_sel_mdl_id, "min.RMSE.OOB"]))

    if (!is.null(glb_category_var)) {
        tmp_OOBobs_df <- glbObsOOB[, c(glb_category_var, glb_rsp_var,
                                           predct_error_var_name)]
        names(tmp_OOBobs_df)[length(names(tmp_OOBobs_df))] <- "error.abs.OOB"
        sOOB_ctgry_df <- dplyr::group_by_(tmp_OOBobs_df, glb_category_var)
        sOOB_ctgry_df <- dplyr::count(sOOB_ctgry_df, 
                                      startprice.log10.abs.OOB.sum = sum(abs(startprice.log10)),
                                        err.abs.OOB.sum = sum(error.abs.OOB),
                                        err.abs.OOB.mean = mean(error.abs.OOB))
        names(sOOB_ctgry_df)[4] <- ".n.OOB"
        sOOB_ctgry_df <- dplyr::ungroup(sOOB_ctgry_df)
        #intersect(names(glb_ctgry_df), names(sOOB_ctgry_df))
        glb_ctgry_df <- merge(glb_ctgry_df, sOOB_ctgry_df, all=TRUE)
        print(orderBy(~-err.abs.OOB.mean, glb_ctgry_df))
    }
    
    if ((glb_rsp_var %in% names(glbObsNew)) &&
        !(any(is.na(glbObsNew[, glb_rsp_var])))) {
            pred_stats_df <- 
                mypredict_mdl(mdl=glb_models_lst[[glb_fin_mdl_id]], 
                              df=glbObsNew, 
                              rsp_var=glb_rsp_var, 
                              rsp_var_out=glb_rsp_var_out, 
                              model_id_method=glb_fin_mdl_id, 
                              label="new",
						      model_summaryFunction=glb_sel_mdl$control$summaryFunction, 
						      model_metric=glb_sel_mdl$metric,
						      model_metric_maximize=glb_sel_mdl$maximize,
						      ret_type="stats")        
            print(sprintf("%s prediction stats for glbObsNew:", glb_fin_mdl_id))
            print(pred_stats_df)
    }    
}    

if (glb_is_classification) {
    print(sprintf("%s OOB confusion matrix & accuracy: ", glb_sel_mdl_id))
    print(t(confusionMatrix(glbObsOOB[, paste0(glb_rsp_var_out, glb_sel_mdl_id)], 
                            glbObsOOB[, glb_rsp_var])$table))

    if (!is.null(glb_category_var)) {
        glb_ctgry_df <- merge(glb_ctgry_df, 
            myget_category_stats(obs_df = glbObsTrn, mdl_id = glb_fin_mdl_id, 
                                 label = "trn"),
                              by = glb_category_var, all = TRUE)
        row.names(glb_ctgry_df) <- glb_ctgry_df[, glb_category_var]
        
        glb_ctgry_df <- merge(glb_ctgry_df, 
            myget_category_stats(obs_df = glbObsNew, mdl_id = glb_fin_mdl_id,
                                 label="new"),
                              by = glb_category_var, all = TRUE)
        row.names(glb_ctgry_df) <- glb_ctgry_df[, glb_category_var]
        
        if (any(grepl("OOB", glbMdlMetricsEval)))
            print(orderBy(~-err.abs.OOB.mean, glb_ctgry_df[, -1])) else
            print(orderBy(~-err.abs.fit.mean, glb_ctgry_df[, -1]))
        print(colSums(glb_ctgry_df[, -grep(glb_category_var,
                                           names(glb_ctgry_df))]))
        
    }
    
    if ((glb_rsp_var %in% names(glbObsNew)) &&
        !(any(is.na(glbObsNew[, glb_rsp_var])))) {
        print(sprintf("%s new confusion matrix & accuracy: ", glb_fin_mdl_id))
        print(t(confusionMatrix(
                            glbObsNew[, paste0(glb_rsp_var_out, glb_fin_mdl_id)], 
                                glbObsNew[, glb_rsp_var])$table))
    }    
}    
```

```
## [1] "All.X.glmnet OOB confusion matrix & accuracy: "
##          Prediction
## Reference    N    Y
##         N 1176  322
##         Y   77  153
##                                    .n.OOB .n.Fit .n.Tst .freqRatio.Fit
## OpEd#Opinion#                          89    437    164    0.090965862
## Styles#U.S.#                           50    127     61    0.026436303
## Science#Health#                        48    148     57    0.030807660
## Business#Crosswords/Games#             18    105     42    0.021856786
## #Opinion#ThePublicEditor                4     16     10    0.003330558
## Business#Technology#                  126    213    114    0.044338052
## Business#BusinessDay#Dealbook         323    629    304    0.130932556
## ##                                    371    913    342    0.190049958
## Metro#N.Y./Region#                     70    128     67    0.026644463
## Culture#Arts#                         185    490    174    0.101998335
## Business#BusinessDay#SmallBusiness     40    100     41    0.020815987
## Styles##Fashion                        15    104     15    0.021648626
## #Opinion#RoomForDebate                 20     42     20    0.008742714
## #Multimedia#                           49     92     52    0.019150708
## Travel#Travel#                         34     83     35    0.017277269
## TStyle##                              101    623    105    0.129683597
## Foreign#World#AsiaPacific              53    150     56    0.031223980
## myOther                                 5     33      5    0.006869276
## Culture##                               1     NA     70             NA
## #U.S.#Education                        82    243     89    0.050582848
## Foreign#World#                         44    128     47    0.026644463
##                                    .freqRatio.OOB .freqRatio.Tst
## OpEd#Opinion#                        0.0515046296    0.087700535
## Styles#U.S.#                         0.0289351852    0.032620321
## Science#Health#                      0.0277777778    0.030481283
## Business#Crosswords/Games#           0.0104166667    0.022459893
## #Opinion#ThePublicEditor             0.0023148148    0.005347594
## Business#Technology#                 0.0729166667    0.060962567
## Business#BusinessDay#Dealbook        0.1869212963    0.162566845
## ##                                   0.2146990741    0.182887701
## Metro#N.Y./Region#                   0.0405092593    0.035828877
## Culture#Arts#                        0.1070601852    0.093048128
## Business#BusinessDay#SmallBusiness   0.0231481481    0.021925134
## Styles##Fashion                      0.0086805556    0.008021390
## #Opinion#RoomForDebate               0.0115740741    0.010695187
## #Multimedia#                         0.0283564815    0.027807487
## Travel#Travel#                       0.0196759259    0.018716578
## TStyle##                             0.0584490741    0.056149733
## Foreign#World#AsiaPacific            0.0306712963    0.029946524
## myOther                              0.0028935185    0.002673797
## Culture##                            0.0005787037    0.037433155
## #U.S.#Education                      0.0474537037    0.047593583
## Foreign#World#                       0.0254629630    0.025133690
##                                    err.abs.fit.sum err.abs.fit.mean .n.fit
## OpEd#Opinion#                           101.962286       0.23332331    437
## Styles#U.S.#                             51.707342       0.40714443    127
## Science#Health#                          48.666469       0.32882749    148
## Business#Crosswords/Games#               19.346059       0.18424818    105
## #Opinion#ThePublicEditor                  4.044658       0.25279116     16
## Business#Technology#                     46.446343       0.21805795    213
## Business#BusinessDay#Dealbook            76.983498       0.12239030    629
## ##                                       91.828632       0.10057901    913
## Metro#N.Y./Region#                       14.444993       0.11285150    128
## Culture#Arts#                            41.599488       0.08489691    490
## Business#BusinessDay#SmallBusiness        8.607239       0.08607239    100
## Styles##Fashion                           3.033172       0.02916512    104
## #Opinion#RoomForDebate                    1.938007       0.04614303     42
## #Multimedia#                              4.739265       0.05151375     92
## Travel#Travel#                            2.569132       0.03095340     83
## TStyle##                                 14.962143       0.02401628    623
## Foreign#World#AsiaPacific                 7.489991       0.04993328    150
## myOther                                   1.720179       0.05212665     33
## Culture##                                       NA               NA     NA
## #U.S.#Education                           4.569704       0.01880537    243
## Foreign#World#                            3.357443       0.02623002    128
##                                    err.abs.OOB.sum err.abs.OOB.mean
## OpEd#Opinion#                          52.11236621       0.58553220
## Styles#U.S.#                           27.39861838       0.54797237
## Science#Health#                        25.12604617       0.52345930
## Business#Crosswords/Games#              9.01410297       0.50078350
## #Opinion#ThePublicEditor                1.97530353       0.49382588
## Business#Technology#                   28.55232120       0.22660572
## Business#BusinessDay#Dealbook          57.75141149       0.17879694
## ##                                     63.39672286       0.17088065
## Metro#N.Y./Region#                     11.12676976       0.15895385
## Culture#Arts#                          28.98681899       0.15668551
## Business#BusinessDay#SmallBusiness      3.84496527       0.09612413
## Styles##Fashion                         1.42711318       0.09514088
## #Opinion#RoomForDebate                  1.78910542       0.08945527
## #Multimedia#                            2.85858104       0.05833839
## Travel#Travel#                          1.97769164       0.05816740
## TStyle##                                5.15834308       0.05107270
## Foreign#World#AsiaPacific               2.69695832       0.05088601
## myOther                                 0.23180569       0.04636114
## Culture##                               0.02849447       0.02849447
## #U.S.#Education                         2.22191005       0.02709646
## Foreign#World#                          1.13916144       0.02589003
##                                    err.abs.trn.sum err.abs.trn.mean .n.trn
## OpEd#Opinion#                         189.20017098       0.35969614    526
## Styles#U.S.#                           83.39117680       0.47113659    177
## Science#Health#                        80.38081529       0.41010620    196
## Business#Crosswords/Games#             36.88161340       0.29985052    123
## #Opinion#ThePublicEditor                7.55559605       0.37777980     20
## Business#Technology#                   76.10301242       0.22449266    339
## Business#BusinessDay#Dealbook         158.25144897       0.16623051    952
## ##                                    180.32449190       0.14043964   1284
## Metro#N.Y./Region#                     30.35960093       0.15333132    198
## Culture#Arts#                          85.20930589       0.12623601    675
## Business#BusinessDay#SmallBusiness     14.98537727       0.10703841    140
## Styles##Fashion                         7.36986388       0.06193163    119
## #Opinion#RoomForDebate                  6.46512194       0.10427616     62
## #Multimedia#                            9.31950181       0.06609576    141
## Travel#Travel#                          6.71907976       0.05742803    117
## TStyle##                               31.23918133       0.04314804    724
## Foreign#World#AsiaPacific              13.72935437       0.06763229    203
## myOther                                 2.82545315       0.07435403     38
## Culture##                               0.05194618       0.05194618      1
## #U.S.#Education                        10.77780084       0.03316246    325
## Foreign#World#                          6.29953592       0.03662521    172
##                                    err.abs.new.sum err.abs.new.mean .n.new
## OpEd#Opinion#                                   NA               NA    164
## Styles#U.S.#                                    NA               NA     61
## Science#Health#                                 NA               NA     57
## Business#Crosswords/Games#                      NA               NA     42
## #Opinion#ThePublicEditor                        NA               NA     10
## Business#Technology#                            NA               NA    114
## Business#BusinessDay#Dealbook                   NA               NA    304
## ##                                              NA               NA    342
## Metro#N.Y./Region#                              NA               NA     67
## Culture#Arts#                                   NA               NA    174
## Business#BusinessDay#SmallBusiness              NA               NA     41
## Styles##Fashion                                 NA               NA     15
## #Opinion#RoomForDebate                          NA               NA     20
## #Multimedia#                                    NA               NA     52
## Travel#Travel#                                  NA               NA     35
## TStyle##                                        NA               NA    105
## Foreign#World#AsiaPacific                       NA               NA     56
## myOther                                         NA               NA      5
## Culture##                                       NA               NA     70
## #U.S.#Education                                 NA               NA     89
## Foreign#World#                                  NA               NA     47
##           .n.OOB           .n.Fit           .n.Tst   .freqRatio.Fit 
##      1728.000000               NA      1870.000000               NA 
##   .freqRatio.OOB   .freqRatio.Tst  err.abs.fit.sum err.abs.fit.mean 
##         1.000000         1.000000               NA               NA 
##           .n.fit  err.abs.OOB.sum err.abs.OOB.mean  err.abs.trn.sum 
##               NA       328.814611         4.170523      1037.439449 
## err.abs.trn.mean           .n.trn  err.abs.new.sum err.abs.new.mean 
##         3.432938      6532.000000               NA               NA 
##           .n.new 
##      1870.000000
```

```r
print(orderBy(as.formula(paste0("~ -", glb_sel_mdl_id, ".importance")), 
              subset(glb_featsimp_df, importance > 10)))
```

```
##                                                    All.X.glmnet.importance
## NDSSName.my.fctrOpEd#Opinion#                                    100.00000
## NDSSName.my.fctrBusiness#Crosswords/Games#                        98.38173
## NDSSName.my.fctr#Opinion#ThePublicEditor                          87.34197
## NDSSName.my.fctrScience#Health#                                   84.05064
## NDSSName.my.fctrStyles#U.S.#                                      80.66804
## NDSSName.my.fctrBusiness#Technology#                              43.29623
## WordCount.log1p                                                   34.01597
## WordCount.root2                                                   32.29795
## .rnorm                                                            31.33944
## NDSSName.my.fctrCulture##                                         31.33944
## NDSSName.my.fctrCulture#Arts#                                     31.33944
## NDSSName.my.fctrMetro#N.Y./Region#                                31.33944
## WordCount.nexp                                                    31.33944
## NDSSName.my.fctrBusiness#BusinessDay#Dealbook                     31.33944
## NDSSName.my.fctrTravel#Travel#                                    31.31570
## NDSSName.my.fctrForeign#World#                                    30.00852
## NDSSName.my.fctr#Multimedia#                                      28.91940
## NDSSName.my.fctrmyOther                                           27.85141
## NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness                27.78546
## NDSSName.my.fctrStyles##Fashion                                   22.02991
## NDSSName.my.fctrForeign#World#AsiaPacific                         19.29994
## NDSSName.my.fctr#U.S.#Education                                   18.25900
## NDSSName.my.fctrTStyle##                                          17.36032
##                                                    importance
## NDSSName.my.fctrOpEd#Opinion#                       100.00000
## NDSSName.my.fctrBusiness#Crosswords/Games#           99.96082
## NDSSName.my.fctr#Opinion#ThePublicEditor             86.61487
## NDSSName.my.fctrScience#Health#                      82.33065
## NDSSName.my.fctrStyles#U.S.#                         77.12194
## NDSSName.my.fctrBusiness#Technology#                 38.62774
## WordCount.log1p                                      37.41260
## WordCount.root2                                      32.88366
## .rnorm                                               31.93272
## NDSSName.my.fctrCulture##                            31.93272
## NDSSName.my.fctrCulture#Arts#                        31.93272
## NDSSName.my.fctrMetro#N.Y./Region#                   31.93272
## WordCount.nexp                                       31.93272
## NDSSName.my.fctrBusiness#BusinessDay#Dealbook        31.12990
## NDSSName.my.fctrTravel#Travel#                       29.62778
## NDSSName.my.fctrForeign#World#                       25.48953
## NDSSName.my.fctr#Multimedia#                         24.73162
## NDSSName.my.fctrmyOther                              26.20874
## NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness   23.59037
## NDSSName.my.fctrStyles##Fashion                      21.10031
## NDSSName.my.fctrForeign#World#AsiaPacific            14.61070
## NDSSName.my.fctr#U.S.#Education                      14.34221
## NDSSName.my.fctrTStyle##                             16.67670
##                                                    Final.glmnet.importance
## NDSSName.my.fctrOpEd#Opinion#                                    100.00000
## NDSSName.my.fctrBusiness#Crosswords/Games#                        99.96082
## NDSSName.my.fctr#Opinion#ThePublicEditor                          86.61487
## NDSSName.my.fctrScience#Health#                                   82.33065
## NDSSName.my.fctrStyles#U.S.#                                      77.12194
## NDSSName.my.fctrBusiness#Technology#                              38.62774
## WordCount.log1p                                                   37.41260
## WordCount.root2                                                   32.88366
## .rnorm                                                            31.93272
## NDSSName.my.fctrCulture##                                         31.93272
## NDSSName.my.fctrCulture#Arts#                                     31.93272
## NDSSName.my.fctrMetro#N.Y./Region#                                31.93272
## WordCount.nexp                                                    31.93272
## NDSSName.my.fctrBusiness#BusinessDay#Dealbook                     31.12990
## NDSSName.my.fctrTravel#Travel#                                    29.62778
## NDSSName.my.fctrForeign#World#                                    25.48953
## NDSSName.my.fctr#Multimedia#                                      24.73162
## NDSSName.my.fctrmyOther                                           26.20874
## NDSSName.my.fctrBusiness#BusinessDay#SmallBusiness                23.59037
## NDSSName.my.fctrStyles##Fashion                                   21.10031
## NDSSName.my.fctrForeign#World#AsiaPacific                         14.61070
## NDSSName.my.fctr#U.S.#Education                                   14.34221
## NDSSName.my.fctrTStyle##                                          16.67670
```

```r
dsp_myCategory_conf_mtrx <- function(myCategory) {
    print(sprintf("%s OOB::myCategory=%s confusion matrix & accuracy: ", 
                  glb_sel_mdl_id, myCategory))
    print(t(confusionMatrix(
        glbObsOOB[glbObsOOB$myCategory == myCategory, 
                      paste0(glb_rsp_var_out, glb_sel_mdl_id)], 
        glbObsOOB[glbObsOOB$myCategory == myCategory, glb_rsp_var])$table))
    print(sum(glbObsOOB[glbObsOOB$myCategory == myCategory, 
                            predct_accurate_var_name]) / 
         nrow(glbObsOOB[glbObsOOB$myCategory == myCategory, ]))
    err_ids <- glbObsOOB[(glbObsOOB$myCategory == myCategory) & 
                             (!glbObsOOB[, predct_accurate_var_name]), glb_id_var]

    OOB_FNerr_df <- glbObsOOB[(glbObsOOB$UniqueID %in% err_ids) & 
                               (glbObsOOB$Popular == 1), 
                        c(
                            ".clusterid", 
                            "Popular", "Headline", "Snippet", "Abstract")]
    print(sprintf("%s OOB::myCategory=%s FN errors: %d", glb_sel_mdl_id, myCategory,
                  nrow(OOB_FNerr_df)))
    print(OOB_FNerr_df)

    OOB_FPerr_df <- glbObsOOB[(glbObsOOB$UniqueID %in% err_ids) & 
                               (glbObsOOB$Popular == 0), 
                        c(
                            ".clusterid", 
                            "Popular", "Headline", "Snippet", "Abstract")]
    print(sprintf("%s OOB::myCategory=%s FP errors: %d", glb_sel_mdl_id, myCategory,
                  nrow(OOB_FPerr_df)))
    print(OOB_FPerr_df)
}
#dsp_myCategory_conf_mtrx(myCategory="OpEd#Opinion#")
#dsp_myCategory_conf_mtrx(myCategory="Business#Business Day#Dealbook")
#dsp_myCategory_conf_mtrx(myCategory="##")

# if (glb_is_classification) {
#     print("FN_OOB_ids:")
#     print(glbObsOOB[glbObsOOB$UniqueID %in% FN_OOB_ids, 
#                         grep(glb_rsp_var, names(glbObsOOB), value=TRUE)])
#     print(glbObsOOB[glbObsOOB$UniqueID %in% FN_OOB_ids, 
#                         glbFeatsText])
#     print(dsp_vctr <- colSums(glbObsOOB[glbObsOOB$UniqueID %in% FN_OOB_ids, 
#                         setdiff(grep("[HSA].", names(glbObsOOB), value=TRUE),
#                                 union(myfind_chr_cols_df(glbObsOOB),
#                     grep(".fctr", names(glbObsOOB), fixed=TRUE, value=TRUE)))]))
# }

print("glbObsNew prediction stats:")
```

```
## [1] "glbObsNew prediction stats:"
```

```r
if (glb_is_regression)
    print(myplot_histogram(glbObsNew, paste0(glb_rsp_var_out, glb_fin_mdl_id)))
if (glb_is_classification)
    print(table(glbObsNew[, paste0(glb_rsp_var_out, glb_fin_mdl_id)]))
```

```
## 
##    N    Y 
## 1144  726
```

```r
# Use this to see how prediction changes by changing one or more values
# players_df <- data.frame(id=c("Chavez", "Giambi", "Menechino", "Myers", "Pena"),
#                          OBP=c(0.338, 0.391, 0.369, 0.313, 0.361),
#                          SLG=c(0.540, 0.450, 0.374, 0.447, 0.500),
#                         cost=c(1400000, 1065000, 295000, 800000, 300000))
# players_df$RS.predict <- predict(glb_models_lst[[csm_mdl_id]], players_df)
# print(orderBy(~ -RS.predict, players_df))
# dsp_chisq.test(Headline.contains="[Vi]deo")

if ((length(diff <- setdiff(names(glbObsTrn), names(glbObsAll))) > 0) ||
    (length(diff <- setdiff(names(glbObsFit), names(glbObsAll))) > 0) ||
    (length(diff <- setdiff(names(glbObsOOB), names(glbObsAll))) > 0) ||
    (length(diff <- setdiff(names(glbObsNew), names(glbObsAll))) > 0)) {
    print(diff)
    stop("glbObs* not in sync")
}

if (glb_save_envir)
    save(glb_feats_df, glbObsAll, 
         #glbObsTrn, glbObsFit, glbObsOOB, glbObsNew,
         glb_models_df, dsp_models_df, glb_models_lst, glb_model_type,
         glb_sel_mdl, glb_sel_mdl_id,
         glb_fin_mdl, glb_fin_mdl_id,
        file=paste0(glb_out_pfx, "prdnew_dsk.RData"))

# tmp_replay_lst <- replay.petrisim(pn=glb_analytics_pn, 
#     replay.trans=(glb_analytics_avl_objs <- c(glb_analytics_avl_objs, 
#         "data.new.prediction")), flip_coord=TRUE)
# print(ggplot.petrinet(tmp_replay_lst[["pn"]]) + coord_flip())

glb_chunks_df <- myadd_chunk(glb_chunks_df, "display.session.info", major.inc=TRUE)
```

```
##                   label step_major step_minor label_minor     bgn     end
## 16     predict.data.new          8          0           0 233.199 242.389
## 17 display.session.info          9          0           0 242.390      NA
##    elapsed
## 16    9.19
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
## 10              fit.models          6          0           0  35.238
## 14       fit.data.training          7          0           0 181.344
## 11              fit.models          6          1           1 139.591
## 12              fit.models          6          2           2 161.330
## 15       fit.data.training          7          1           1 222.039
## 9          select.features          5          0           0  25.428
## 1              import.data          1          0           0   9.037
## 16        predict.data.new          8          0           0 233.199
## 13              fit.models          6          3           3 174.894
## 2             inspect.data          2          0           0  18.458
## 5         extract.features          3          0           0  22.532
## 8  partition.data.training          4          0           0  24.288
## 3               scrub.data          2          1           1  21.340
## 4           transform.data          2          2           2  22.441
## 6      manage.missing.data          3          1           1  24.178
## 7             cluster.data          3          2           2  24.259
##        end elapsed duration
## 10 139.590 104.352  104.352
## 14 222.039  40.695   40.695
## 11 161.329  21.738   21.738
## 12 174.893  13.563   13.563
## 15 233.199  11.160   11.160
## 9   35.237   9.809    9.809
## 1   18.458   9.421    9.421
## 16 242.389   9.190    9.190
## 13 181.343   6.450    6.449
## 2   21.340   2.882    2.882
## 5   24.177   1.645    1.645
## 8   25.428   1.140    1.140
## 3   22.440   1.100    1.100
## 4   22.531   0.091    0.090
## 6   24.259   0.081    0.081
## 7   24.288   0.029    0.029
## [1] "Total Elapsed Time: 242.389 secs"
```

![](NYTBlogs3_base_files/figure-html/display.session.info-1.png) 
