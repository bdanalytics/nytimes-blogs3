# NYTimes:Blogs: Popular classification:: NYTBlogs3_feat_PubDate
bdanalytics  

**  **    
**Date: (Sat) Nov 21, 2015**    

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
#glbFeatsInteractionOnly[["carrier.fctr"]] <- "cellular.fctr"

# currently does not handle more than 1 column; consider concatenating multiple columns
glb_id_var <- "UniqueID" # choose from c(NULL : default, "<id_feat>") 
glbFeatsCategory <- "NDSSName.my.fctr" # choose from c(NULL : default, "<category>")

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
glbFeatsDerive[["NDSSName.my"]] <- list(
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
    c(format = "%Y-%m-%d %H:%M:%S", timezone = "America/New_York", impute.na = TRUE, 
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
glb_tune_models_df <- data.frame()
# Experiment specific code to avoid caret crash
glmnet_tune_models_df <- rbind(data.frame()
                            ,data.frame(method = "glmnet", parameter = "alpha", 
                                        vals = "0.100 0.325 0.550 0.775 1.000")
                            ,data.frame(method = "glmnet", parameter = "lambda",
                                        vals = "9.342e-02")    
                                    )

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
##         label step_major step_minor label_minor   bgn end elapsed
## 1 import.data          1          0           0 9.172  NA      NA
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

![](NYTBlogs3_feat_PubDate_files/figure-html/import.data-2.png) 

```
##    .src   .n
## 1 Train 6532
## 2  Test 1870
```

```
##          label step_major step_minor label_minor    bgn    end elapsed
## 1  import.data          1          0           0  9.172 20.291  11.119
## 2 inspect.data          2          0           0 20.292     NA      NA
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
```

![](NYTBlogs3_feat_PubDate_files/figure-html/inspect.data-3.png) 

```
##          label step_major step_minor label_minor    bgn    end elapsed
## 2 inspect.data          2          0           0 20.292 23.326   3.034
## 3   scrub.data          2          1           1 23.327     NA      NA
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
## 3     scrub.data          2          1           1 23.327 24.356   1.029
## 4 transform.data          2          2           2 24.356     NA      NA
```

### Step `2.2: transform data`

```
## [1] "Creating new feature: NDSSName.my..."
## [1] "Creating new feature: WordCount.log1p..."
## [1] "Creating new feature: WordCount.root2..."
## [1] "Creating new feature: WordCount.nexp..."
```

```
##              label step_major step_minor label_minor    bgn    end elapsed
## 4   transform.data          2          2           2 24.356 24.452   0.096
## 5 extract.features          3          0           0 24.453     NA      NA
```

## Step `3.0: extract features`

```
##                  label step_major step_minor label_minor    bgn end
## 1 extract.features_bgn          1          0           0 24.509  NA
##   elapsed
## 1      NA
```

```
##                                 label step_major step_minor label_minor
## 1                extract.features_bgn          1          0           0
## 2 extract.features_factorize.str.vars          2          0           0
##      bgn    end elapsed
## 1 24.509 24.519    0.01
## 2 24.520     NA      NA
```

```
##         NewsDesk      SectionName   SubsectionName         Headline 
##       "NewsDesk"    "SectionName" "SubsectionName"       "Headline" 
##          Snippet         Abstract          PubDate             .src 
##        "Snippet"       "Abstract"        "PubDate"           ".src" 
##      NDSSName.my 
##    "NDSSName.my"
```

```
## Warning: Creating factors of string variable: NDSSName.my: # of unique
## values: 21
```

```
##                                   label step_major step_minor label_minor
## 2   extract.features_factorize.str.vars          2          0           0
## 3 extract.features_xtract.DateTime.vars          3          0           0
##      bgn    end elapsed
## 2 24.520 24.537   0.017
## 3 24.538     NA      NA
## [1] "Extracting features from DateTime(s): PubDate"
```

```
## Loading required package: XML
```

```
## [1] "**********"
## [1] "Consider adding state & city holidays for glbFeatsDateTime: PubDate"
## [1] "**********"
```

```
## Loading required package: zoo
## 
## Attaching package: 'zoo'
## 
## The following objects are masked from 'package:base':
## 
##     as.Date, as.Date.numeric
```

```
## [1] "Missing data for numerics:"
##  PubDate.last2.log1p  PubDate.last4.log1p  PubDate.last8.log1p 
##                    2                    4                    8 
## PubDate.last16.log1p PubDate.last32.log1p 
##                   16                   32
```

```
## Loading required package: mice
## Loading required package: Rcpp
## mice 2.25 2015-11-09
```

![](NYTBlogs3_feat_PubDate_files/figure-html/extract.features-1.png) 

```
## [1] "Summary before imputation: "
##    NewsDesk         SectionName        SubsectionName    
##  Length:8402        Length:8402        Length:8402       
##  Class :character   Class :character   Class :character  
##  Mode  :character   Mode  :character   Mode  :character  
##                                                          
##                                                          
##                                                          
##                                                          
##    Headline           Snippet            Abstract        
##  Length:8402        Length:8402        Length:8402       
##  Class :character   Class :character   Class :character  
##  Mode  :character   Mode  :character   Mode  :character  
##                                                          
##                                                          
##                                                          
##                                                          
##    WordCount         PubDate             UniqueID        .src          
##  Min.   :    0.0   Length:8402        Min.   :   1   Length:8402       
##  1st Qu.:  188.0   Class :character   1st Qu.:2101   Class :character  
##  Median :  377.0   Mode  :character   Median :4202   Mode  :character  
##  Mean   :  528.8                      Mean   :4202                     
##  3rd Qu.:  735.0                      3rd Qu.:6302                     
##  Max.   :10912.0                      Max.   :8402                     
##                                                                        
##      .rnorm         NDSSName.my        WordCount.log1p WordCount.root2 
##  Min.   :-3.89398   Length:8402        Min.   :0.000   Min.   :  0.00  
##  1st Qu.:-0.65896   Class :character   1st Qu.:5.242   1st Qu.: 13.71  
##  Median : 0.02058   Mode  :character   Median :5.935   Median : 19.42  
##  Mean   : 0.01549                      Mean   :5.751   Mean   : 20.60  
##  3rd Qu.: 0.67764                      3rd Qu.:6.601   3rd Qu.: 27.11  
##  Max.   : 3.74468                      Max.   :9.298   Max.   :104.46  
##                                                                        
##  WordCount.nexp              NDSSName.my.fctr PubDate.year.fctr
##  Min.   :0.00000   ##                :1626    2014:8402        
##  1st Qu.:0.00000   Bsnss#BsnssDy#Dlbk:1256                     
##  Median :0.00000   Cltr#Arts#        : 849                     
##  Mean   :0.01318   TStyl##           : 829                     
##  3rd Qu.:0.00000   OpEd#Opnn#        : 690                     
##  Max.   :1.00000   Bsnss#Tchnlgy#    : 453                     
##                    (Other)           :2699                     
##  PubDate.month.fctr PubDate.date.fctr PubDate.juliandate
##  09:2341            (0.97,7]:1981     Min.   :244.0     
##  10:2382            (7,13]  :1757     1st Qu.:270.0     
##  11:1809            (13,19] :1808     Median :297.0     
##  12:1870            (19,25] :1650     Mean   :300.1     
##                     (25,31] :1206     3rd Qu.:328.0     
##                                       Max.   :365.0     
##                                                         
##  PubDate.wkday.fctr PubDate.wkend    PubDate.hlday   
##  0: 378             Min.   :0.0000   Min.   :0.0000  
##  1:1605             1st Qu.:0.0000   1st Qu.:0.0000  
##  2:1559             Median :0.0000   Median :0.0000  
##  3:1614             Mean   :0.0732   Mean   :0.0288  
##  4:1539             3rd Qu.:0.0000   3rd Qu.:0.0000  
##  5:1470             Max.   :1.0000   Max.   :1.0000  
##  6: 237                                              
##      PubDate.hour.fctr    PubDate.minute.fctr    PubDate.second.fctr
##  (-0.023,7.67]:1610    (-0.059,14.8]:3119     (-0.059,14.8]:2134    
##  (7.67,15.3]  :4484    (14.8,29.5]  :1671     (14.8,29.5]  :2063    
##  (15.3,23]    :2308    (29.5,44.2]  :1995     (29.5,44.2]  :2112    
##                        (44.2,59.1]  :1617     (44.2,59.1]  :2093    
##                                                                     
##                                                                     
##                                                                     
##  PubDate.day.minutes PubDate.day.minutes.poly.1 PubDate.day.minutes.poly.2
##  Min.   :   0.0      Min.   :-0.0274946         Min.   :-0.008759         
##  1st Qu.: 540.0      1st Qu.:-0.0078858         1st Qu.:-0.007794         
##  Median : 765.0      Median : 0.0002845         Median :-0.003965         
##  Mean   : 757.2      Mean   : 0.0000000         Mean   : 0.000000         
##  3rd Qu.: 977.8      3rd Qu.: 0.0080100         3rd Qu.: 0.002979         
##  Max.   :1439.0      Max.   : 0.0247592         Max.   : 0.042684         
##                                                                           
##  PubDate.day.minutes.poly.3 PubDate.day.minutes.poly.4
##  Min.   :-0.045125          Min.   :-0.0183274        
##  1st Qu.:-0.007780          1st Qu.:-0.0085001        
##  Median :-0.000337          Median : 0.0007244        
##  Mean   : 0.000000          Mean   : 0.0000000        
##  3rd Qu.: 0.007554          3rd Qu.: 0.0081632        
##  Max.   : 0.052153          Max.   : 0.0667744        
##                                                       
##  PubDate.day.minutes.poly.5     .order     PubDate.last2.log1p
##  Min.   :-0.0245092         Min.   :   1   Min.   : 0.6931    
##  1st Qu.:-0.0082017         1st Qu.:2101   1st Qu.: 6.5392    
##  Median :-0.0003821         Median :4202   Median : 7.1835    
##  Mean   : 0.0000000         Mean   :4202   Mean   : 7.1698    
##  3rd Qu.: 0.0077006         3rd Qu.:6302   3rd Qu.: 7.7932    
##  Max.   : 0.0847176         Max.   :8402   Max.   :10.9120    
##                                            NA's   :2          
##  PubDate.last4.log1p PubDate.last8.log1p PubDate.last16.log1p
##  Min.   : 4.159      Min.   : 6.667      Min.   : 8.001      
##  1st Qu.: 7.484      1st Qu.: 8.259      1st Qu.: 9.030      
##  Median : 7.923      Median : 8.626      Median : 9.343      
##  Mean   : 8.039      Mean   : 8.817      Mean   : 9.567      
##  3rd Qu.: 8.465      3rd Qu.: 9.191      3rd Qu.:10.126      
##  Max.   :11.243      Max.   :11.622      Max.   :11.957      
##  NA's   :4           NA's   :8           NA's   :16          
##  PubDate.last32.log1p
##  Min.   : 8.836      
##  1st Qu.: 9.784      
##  Median :10.121      
##  Mean   :10.331      
##  3rd Qu.:10.801      
##  Max.   :12.323      
##  NA's   :32          
## 
##  iter imp variable
##   1   1  PubDate.last2.log1p  PubDate.last4.log1p  PubDate.last8.log1p  PubDate.last16.log1p  PubDate.last32.log1p
##   1   2  PubDate.last2.log1p  PubDate.last4.log1p  PubDate.last8.log1p  PubDate.last16.log1p  PubDate.last32.log1p
##   1   3  PubDate.last2.log1p  PubDate.last4.log1p  PubDate.last8.log1p  PubDate.last16.log1p  PubDate.last32.log1p
##   1   4  PubDate.last2.log1p  PubDate.last4.log1p  PubDate.last8.log1p  PubDate.last16.log1p  PubDate.last32.log1p
##   1   5  PubDate.last2.log1p  PubDate.last4.log1p  PubDate.last8.log1p  PubDate.last16.log1p  PubDate.last32.log1p
##   2   1  PubDate.last2.log1p  PubDate.last4.log1p  PubDate.last8.log1p  PubDate.last16.log1p  PubDate.last32.log1p
##   2   2  PubDate.last2.log1p  PubDate.last4.log1p  PubDate.last8.log1p  PubDate.last16.log1p  PubDate.last32.log1p
##   2   3  PubDate.last2.log1p  PubDate.last4.log1p  PubDate.last8.log1p  PubDate.last16.log1p  PubDate.last32.log1p
##   2   4  PubDate.last2.log1p  PubDate.last4.log1p  PubDate.last8.log1p  PubDate.last16.log1p  PubDate.last32.log1p
##   2   5  PubDate.last2.log1p  PubDate.last4.log1p  PubDate.last8.log1p  PubDate.last16.log1p  PubDate.last32.log1p
##   3   1  PubDate.last2.log1p  PubDate.last4.log1p  PubDate.last8.log1p  PubDate.last16.log1p  PubDate.last32.log1p
##   3   2  PubDate.last2.log1p  PubDate.last4.log1p  PubDate.last8.log1p  PubDate.last16.log1p  PubDate.last32.log1p
##   3   3  PubDate.last2.log1p  PubDate.last4.log1p  PubDate.last8.log1p  PubDate.last16.log1p  PubDate.last32.log1p
##   3   4  PubDate.last2.log1p  PubDate.last4.log1p  PubDate.last8.log1p  PubDate.last16.log1p  PubDate.last32.log1p
##   3   5  PubDate.last2.log1p  PubDate.last4.log1p  PubDate.last8.log1p  PubDate.last16.log1p  PubDate.last32.log1p
##   4   1  PubDate.last2.log1p  PubDate.last4.log1p  PubDate.last8.log1p  PubDate.last16.log1p  PubDate.last32.log1p
##   4   2  PubDate.last2.log1p  PubDate.last4.log1p  PubDate.last8.log1p  PubDate.last16.log1p  PubDate.last32.log1p
##   4   3  PubDate.last2.log1p  PubDate.last4.log1p  PubDate.last8.log1p  PubDate.last16.log1p  PubDate.last32.log1p
##   4   4  PubDate.last2.log1p  PubDate.last4.log1p  PubDate.last8.log1p  PubDate.last16.log1p  PubDate.last32.log1p
##   4   5  PubDate.last2.log1p  PubDate.last4.log1p  PubDate.last8.log1p  PubDate.last16.log1p  PubDate.last32.log1p
##   5   1  PubDate.last2.log1p  PubDate.last4.log1p  PubDate.last8.log1p  PubDate.last16.log1p  PubDate.last32.log1p
##   5   2  PubDate.last2.log1p  PubDate.last4.log1p  PubDate.last8.log1p  PubDate.last16.log1p  PubDate.last32.log1p
##   5   3  PubDate.last2.log1p  PubDate.last4.log1p  PubDate.last8.log1p  PubDate.last16.log1p  PubDate.last32.log1p
##   5   4  PubDate.last2.log1p  PubDate.last4.log1p  PubDate.last8.log1p  PubDate.last16.log1p  PubDate.last32.log1p
##   5   5  PubDate.last2.log1p  PubDate.last4.log1p  PubDate.last8.log1p  PubDate.last16.log1p  PubDate.last32.log1p
##    WordCount          UniqueID        .rnorm         WordCount.log1p
##  Min.   :    0.0   Min.   :   1   Min.   :-3.89398   Min.   :0.000  
##  1st Qu.:  188.0   1st Qu.:2101   1st Qu.:-0.65896   1st Qu.:5.242  
##  Median :  377.0   Median :4202   Median : 0.02058   Median :5.935  
##  Mean   :  528.8   Mean   :4202   Mean   : 0.01549   Mean   :5.751  
##  3rd Qu.:  735.0   3rd Qu.:6302   3rd Qu.: 0.67764   3rd Qu.:6.601  
##  Max.   :10912.0   Max.   :8402   Max.   : 3.74468   Max.   :9.298  
##                                                                     
##  WordCount.root2  WordCount.nexp              NDSSName.my.fctr
##  Min.   :  0.00   Min.   :0.00000   ##                :1626   
##  1st Qu.: 13.71   1st Qu.:0.00000   Bsnss#BsnssDy#Dlbk:1256   
##  Median : 19.42   Median :0.00000   Cltr#Arts#        : 849   
##  Mean   : 20.60   Mean   :0.01318   TStyl##           : 829   
##  3rd Qu.: 27.11   3rd Qu.:0.00000   OpEd#Opnn#        : 690   
##  Max.   :104.46   Max.   :1.00000   Bsnss#Tchnlgy#    : 453   
##                                     (Other)           :2699   
##  PubDate.year.fctr PubDate.month.fctr PubDate.date.fctr PubDate.juliandate
##  2014:8402         09:2341            (0.97,7]:1981     Min.   :244.0     
##                    10:2382            (7,13]  :1757     1st Qu.:270.0     
##                    11:1809            (13,19] :1808     Median :297.0     
##                    12:1870            (19,25] :1650     Mean   :300.1     
##                                       (25,31] :1206     3rd Qu.:328.0     
##                                                         Max.   :365.0     
##                                                                           
##  PubDate.wkday.fctr PubDate.wkend    PubDate.hlday   
##  0: 378             Min.   :0.0000   Min.   :0.0000  
##  1:1605             1st Qu.:0.0000   1st Qu.:0.0000  
##  2:1559             Median :0.0000   Median :0.0000  
##  3:1614             Mean   :0.0732   Mean   :0.0288  
##  4:1539             3rd Qu.:0.0000   3rd Qu.:0.0000  
##  5:1470             Max.   :1.0000   Max.   :1.0000  
##  6: 237                                              
##      PubDate.hour.fctr    PubDate.minute.fctr    PubDate.second.fctr
##  (-0.023,7.67]:1610    (-0.059,14.8]:3119     (-0.059,14.8]:2134    
##  (7.67,15.3]  :4484    (14.8,29.5]  :1671     (14.8,29.5]  :2063    
##  (15.3,23]    :2308    (29.5,44.2]  :1995     (29.5,44.2]  :2112    
##                        (44.2,59.1]  :1617     (44.2,59.1]  :2093    
##                                                                     
##                                                                     
##                                                                     
##  PubDate.day.minutes PubDate.day.minutes.poly.1 PubDate.day.minutes.poly.2
##  Min.   :   0.0      Min.   :-0.0274946         Min.   :-0.008759         
##  1st Qu.: 540.0      1st Qu.:-0.0078858         1st Qu.:-0.007794         
##  Median : 765.0      Median : 0.0002845         Median :-0.003965         
##  Mean   : 757.2      Mean   : 0.0000000         Mean   : 0.000000         
##  3rd Qu.: 977.8      3rd Qu.: 0.0080100         3rd Qu.: 0.002979         
##  Max.   :1439.0      Max.   : 0.0247592         Max.   : 0.042684         
##                                                                           
##  PubDate.day.minutes.poly.3 PubDate.day.minutes.poly.4
##  Min.   :-0.045125          Min.   :-0.0183274        
##  1st Qu.:-0.007780          1st Qu.:-0.0085001        
##  Median :-0.000337          Median : 0.0007244        
##  Mean   : 0.000000          Mean   : 0.0000000        
##  3rd Qu.: 0.007554          3rd Qu.: 0.0081632        
##  Max.   : 0.052153          Max.   : 0.0667744        
##                                                       
##  PubDate.day.minutes.poly.5     .order     PubDate.last2.log1p
##  Min.   :-0.0245092         Min.   :   1   Min.   : 0.6931    
##  1st Qu.:-0.0082017         1st Qu.:2101   1st Qu.: 6.5396    
##  Median :-0.0003821         Median :4202   Median : 7.1839    
##  Mean   : 0.0000000         Mean   :4202   Mean   : 7.1701    
##  3rd Qu.: 0.0077006         3rd Qu.:6302   3rd Qu.: 7.7941    
##  Max.   : 0.0847176         Max.   :8402   Max.   :10.9120    
##                                                               
##  PubDate.last4.log1p PubDate.last8.log1p PubDate.last16.log1p
##  Min.   : 4.159      Min.   : 6.667      Min.   : 8.001      
##  1st Qu.: 7.484      1st Qu.: 8.259      1st Qu.: 9.030      
##  Median : 7.923      Median : 8.627      Median : 9.344      
##  Mean   : 8.040      Mean   : 8.818      Mean   : 9.569      
##  3rd Qu.: 8.466      3rd Qu.: 9.193      3rd Qu.:10.130      
##  Max.   :11.243      Max.   :11.622      Max.   :11.957      
##                                                              
##  PubDate.last32.log1p
##  Min.   : 8.836      
##  1st Qu.: 9.785      
##  Median :10.125      
##  Mean   :10.334      
##  3rd Qu.:10.804      
##  Max.   :12.323      
## 
```

```
##                                   label step_major step_minor label_minor
## 3 extract.features_xtract.DateTime.vars          3          0           0
## 4                  extract.features_end          4          0           0
##      bgn    end elapsed
## 3 24.538 58.264  33.726
## 4 58.265     NA      NA
```

```
##                                   label step_major step_minor label_minor
## 3 extract.features_xtract.DateTime.vars          3          0           0
## 2   extract.features_factorize.str.vars          2          0           0
## 1                  extract.features_bgn          1          0           0
##      bgn    end elapsed duration
## 3 24.538 58.264  33.726   33.726
## 2 24.520 24.537   0.017    0.017
## 1 24.509 24.519   0.010    0.010
## [1] "Total Elapsed Time: 58.264 secs"
```

![](NYTBlogs3_feat_PubDate_files/figure-html/extract.features-2.png) 

```
## time	trans	 "bgn " "fit.data.training.all " "predict.data.new " "end " 
## 0.0000 	multiple enabled transitions:  data.training.all data.new model.selected 	firing:  data.training.all 
## 1.0000 	 1 	 2 1 0 0 
## 1.0000 	multiple enabled transitions:  data.training.all data.new model.selected model.final data.training.all.prediction 	firing:  data.new 
## 2.0000 	 2 	 1 1 1 0
```

![](NYTBlogs3_feat_PubDate_files/figure-html/extract.features-3.png) 

```
##                 label step_major step_minor label_minor    bgn    end
## 5    extract.features          3          0           0 24.453 59.576
## 6 manage.missing.data          3          1           1 59.576     NA
##   elapsed
## 5  35.123
## 6      NA
```

### Step `3.1: manage missing data`

```
## [1] "numeric data missing in : "
##   Popular Pplr.fctr 
##      1870      1870 
## [1] "numeric data w/ 0s in : "
##           WordCount             Popular     WordCount.log1p 
##                 109                5439                 109 
##     WordCount.root2      WordCount.nexp  PubDate.wkday.fctr 
##                 109                2044                 378 
##       PubDate.wkend       PubDate.hlday PubDate.day.minutes 
##                7787                8160                   5 
## [1] "numeric data w/ Infs in : "
## named integer(0)
## [1] "numeric data w/ NaNs in : "
## named integer(0)
## [1] "string data missing in : "
##       NewsDesk    SectionName SubsectionName       Headline        Snippet 
##           2408           2899           6176              0             13 
##       Abstract        PubDate    NDSSName.my 
##             17              0              0
```

```
## [1] "numeric data missing in : "
##   Popular Pplr.fctr 
##      1870      1870 
## [1] "numeric data w/ 0s in : "
##           WordCount             Popular     WordCount.log1p 
##                 109                5439                 109 
##     WordCount.root2      WordCount.nexp  PubDate.wkday.fctr 
##                 109                2044                 378 
##       PubDate.wkend       PubDate.hlday PubDate.day.minutes 
##                7787                8160                   5 
## [1] "numeric data w/ Infs in : "
## named integer(0)
## [1] "numeric data w/ NaNs in : "
## named integer(0)
## [1] "string data missing in : "
##       NewsDesk    SectionName SubsectionName       Headline        Snippet 
##           2408           2899           6176              0             13 
##       Abstract        PubDate    NDSSName.my 
##             17              0              0
```

```
##                 label step_major step_minor label_minor    bgn    end
## 6 manage.missing.data          3          1           1 59.576 60.676
## 7        cluster.data          3          2           2 60.677     NA
##   elapsed
## 6   1.101
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
## 7            cluster.data          3          2           2 60.677 60.745
## 8 partition.data.training          4          0           0 60.746     NA
##   elapsed
## 7   0.068
## 8      NA
```

## Step `4.0: partition data training`

```
## [1] "Prediction Hints by Catgeory:"
##    NDSSName.my.fctr Popular.0 Popular.1 .n.tst .strata.0 .strata.1
## 5       #U.S.#Edctn       325        NA     89        82        17
## 10           Cltr##         1        NA     70         1        13
## 12       Frgn#Wrld#       172        NA     47        44         9
## 21           myOthr        38        NA      5         5         1
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
##           NDSSName.my.fctr .n.Fit .n.OOB .n.Tst .freqRatio.Fit
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
##                     label step_major step_minor label_minor    bgn   end
## 8 partition.data.training          4          0           0 60.746 62.06
## 9         select.features          5          0           0 62.061    NA
##   elapsed
## 8   1.314
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
## NDSSName.my.fctr                     NDSSName.my.fctr  0.165445970
## PubDate.day.minutes               PubDate.day.minutes  0.156753478
## PubDate.day.minutes.poly.1 PubDate.day.minutes.poly.1  0.156753478
## PubDate.hour.fctr                   PubDate.hour.fctr  0.135436805
## PubDate.wkend                           PubDate.wkend  0.104707290
## PubDate.day.minutes.poly.4 PubDate.day.minutes.poly.4  0.073941394
## PubDate.day.minutes.poly.2 PubDate.day.minutes.poly.2  0.070977720
## PubDate.last4.log1p               PubDate.last4.log1p  0.069776398
## PubDate.last2.log1p               PubDate.last2.log1p  0.065679536
## PubDate.last8.log1p               PubDate.last8.log1p  0.056574578
## PubDate.day.minutes.poly.5 PubDate.day.minutes.poly.5 -0.055929231
## WordCount.nexp                         WordCount.nexp -0.053208396
## PubDate.wkday.fctr                 PubDate.wkday.fctr -0.039801288
## PubDate.last16.log1p             PubDate.last16.log1p  0.038456811
## PubDate.minute.fctr               PubDate.minute.fctr -0.034073846
## PubDate.day.minutes.poly.3 PubDate.day.minutes.poly.3  0.027983551
## PubDate.last32.log1p             PubDate.last32.log1p  0.022542505
## PubDate.month.fctr                 PubDate.month.fctr  0.019148739
## PubDate.POSIX                           PubDate.POSIX  0.015683258
## PubDate.hlday                           PubDate.hlday  0.014690122
## PubDate.juliandate                 PubDate.juliandate  0.014361075
## PubDate.zoo                               PubDate.zoo  0.013260902
## PubDate.second.fctr               PubDate.second.fctr -0.011879458
## UniqueID                                     UniqueID  0.011824920
## PubDate.date.fctr                   PubDate.date.fctr -0.011647558
## .rnorm                                         .rnorm  0.008212201
## PubDate.year.fctr                   PubDate.year.fctr           NA
##                            exclude.as.feat   cor.y.abs
## Popular                                  1 1.000000000
## WordCount.root2                          0 0.292120679
## WordCount                                1 0.257526549
## WordCount.log1p                          0 0.254319628
## NDSSName.my.fctr                         0 0.165445970
## PubDate.day.minutes                      1 0.156753478
## PubDate.day.minutes.poly.1               0 0.156753478
## PubDate.hour.fctr                        0 0.135436805
## PubDate.wkend                            0 0.104707290
## PubDate.day.minutes.poly.4               0 0.073941394
## PubDate.day.minutes.poly.2               0 0.070977720
## PubDate.last4.log1p                      0 0.069776398
## PubDate.last2.log1p                      0 0.065679536
## PubDate.last8.log1p                      0 0.056574578
## PubDate.day.minutes.poly.5               0 0.055929231
## WordCount.nexp                           0 0.053208396
## PubDate.wkday.fctr                       0 0.039801288
## PubDate.last16.log1p                     0 0.038456811
## PubDate.minute.fctr                      0 0.034073846
## PubDate.day.minutes.poly.3               0 0.027983551
## PubDate.last32.log1p                     0 0.022542505
## PubDate.month.fctr                       0 0.019148739
## PubDate.POSIX                            1 0.015683258
## PubDate.hlday                            0 0.014690122
## PubDate.juliandate                       0 0.014361075
## PubDate.zoo                              1 0.013260902
## PubDate.second.fctr                      0 0.011879458
## UniqueID                                 1 0.011824920
## PubDate.date.fctr                        0 0.011647558
## .rnorm                                   0 0.008212201
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
## [1] "cor(PubDate.last16.log1p, PubDate.last8.log1p)=0.8942"
## [1] "cor(Pplr.fctr, PubDate.last16.log1p)=0.0385"
## [1] "cor(Pplr.fctr, PubDate.last8.log1p)=0.0566"
```

```
## Warning in myfind_cor_features(feats_df = glb_feats_df, obs_df =
## glbObsTrn, : Identified PubDate.last16.log1p as highly correlated with
## PubDate.last8.log1p
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
## [1] "cor(PubDate.last4.log1p, PubDate.last8.log1p)=0.8579"
## [1] "cor(Pplr.fctr, PubDate.last4.log1p)=0.0698"
## [1] "cor(Pplr.fctr, PubDate.last8.log1p)=0.0566"
```

```
## Warning in myfind_cor_features(feats_df = glb_feats_df, obs_df =
## glbObsTrn, : Identified PubDate.last8.log1p as highly correlated with
## PubDate.last4.log1p
```

```
## [1] "cor(PubDate.last2.log1p, PubDate.last4.log1p)=0.7709"
## [1] "cor(Pplr.fctr, PubDate.last2.log1p)=0.0657"
## [1] "cor(Pplr.fctr, PubDate.last4.log1p)=0.0698"
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
## NDSSName.my.fctr                     NDSSName.my.fctr  0.165445970
## PubDate.day.minutes               PubDate.day.minutes  0.156753478
## PubDate.day.minutes.poly.1 PubDate.day.minutes.poly.1  0.156753478
## PubDate.hour.fctr                   PubDate.hour.fctr  0.135436805
## PubDate.wkend                           PubDate.wkend  0.104707290
## PubDate.day.minutes.poly.4 PubDate.day.minutes.poly.4  0.073941394
## PubDate.day.minutes.poly.2 PubDate.day.minutes.poly.2  0.070977720
## PubDate.last4.log1p               PubDate.last4.log1p  0.069776398
## PubDate.last2.log1p               PubDate.last2.log1p  0.065679536
## PubDate.last8.log1p               PubDate.last8.log1p  0.056574578
## PubDate.last16.log1p             PubDate.last16.log1p  0.038456811
## PubDate.day.minutes.poly.3 PubDate.day.minutes.poly.3  0.027983551
## PubDate.last32.log1p             PubDate.last32.log1p  0.022542505
## PubDate.month.fctr                 PubDate.month.fctr  0.019148739
## PubDate.POSIX                           PubDate.POSIX  0.015683258
## PubDate.hlday                           PubDate.hlday  0.014690122
## PubDate.juliandate                 PubDate.juliandate  0.014361075
## PubDate.zoo                               PubDate.zoo  0.013260902
## UniqueID                                     UniqueID  0.011824920
## .rnorm                                         .rnorm  0.008212201
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
## NDSSName.my.fctr                         0 0.165445970
## PubDate.day.minutes                      1 0.156753478
## PubDate.day.minutes.poly.1               0 0.156753478
## PubDate.hour.fctr                        0 0.135436805
## PubDate.wkend                            0 0.104707290
## PubDate.day.minutes.poly.4               0 0.073941394
## PubDate.day.minutes.poly.2               0 0.070977720
## PubDate.last4.log1p                      0 0.069776398
## PubDate.last2.log1p                      0 0.065679536
## PubDate.last8.log1p                      0 0.056574578
## PubDate.last16.log1p                     0 0.038456811
## PubDate.day.minutes.poly.3               0 0.027983551
## PubDate.last32.log1p                     0 0.022542505
## PubDate.month.fctr                       0 0.019148739
## PubDate.POSIX                            1 0.015683258
## PubDate.hlday                            0 0.014690122
## PubDate.juliandate                       0 0.014361075
## PubDate.zoo                              1 0.013260902
## UniqueID                                 1 0.011824920
## .rnorm                                   0 0.008212201
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
## NDSSName.my.fctr                                 <NA>  1.348739
## PubDate.day.minutes                              <NA>  1.225490
## PubDate.day.minutes.poly.1                       <NA>  1.225490
## PubDate.hour.fctr          PubDate.day.minutes.poly.1  1.835040
## PubDate.wkend                                    <NA> 12.011952
## PubDate.day.minutes.poly.4                       <NA>  1.225490
## PubDate.day.minutes.poly.2                       <NA>  1.225490
## PubDate.last4.log1p                              <NA>  1.125000
## PubDate.last2.log1p               PubDate.last4.log1p  1.375000
## PubDate.last8.log1p               PubDate.last4.log1p  1.166667
## PubDate.last16.log1p              PubDate.last8.log1p  1.000000
## PubDate.day.minutes.poly.3                       <NA>  1.225490
## PubDate.last32.log1p                             <NA>  1.000000
## PubDate.month.fctr                               <NA>  1.017514
## PubDate.POSIX                                    <NA>  1.000000
## PubDate.hlday                                    <NA> 28.160714
## PubDate.juliandate                 PubDate.month.fctr  1.032520
## PubDate.zoo                                      <NA>  1.000000
## UniqueID                                         <NA>  1.000000
## .rnorm                                           <NA>  1.000000
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
## NDSSName.my.fctr              0.32149418   FALSE FALSE            FALSE
## PubDate.day.minutes          18.08022045   FALSE FALSE            FALSE
## PubDate.day.minutes.poly.1   18.08022045   FALSE FALSE            FALSE
## PubDate.hour.fctr             0.04592774   FALSE FALSE            FALSE
## PubDate.wkend                 0.03061849   FALSE FALSE            FALSE
## PubDate.day.minutes.poly.4   18.08022045   FALSE FALSE            FALSE
## PubDate.day.minutes.poly.2   18.08022045   FALSE FALSE            FALSE
## PubDate.last4.log1p          64.98775260   FALSE FALSE            FALSE
## PubDate.last2.log1p          51.16350276   FALSE FALSE            FALSE
## PubDate.last8.log1p          75.15309247   FALSE FALSE            FALSE
## PubDate.last16.log1p         84.50704225   FALSE FALSE            FALSE
## PubDate.day.minutes.poly.3   18.08022045   FALSE FALSE            FALSE
## PubDate.last32.log1p         91.13594611   FALSE FALSE            FALSE
## PubDate.month.fctr            0.04592774   FALSE FALSE            FALSE
## PubDate.POSIX                99.86221678   FALSE FALSE            FALSE
## PubDate.hlday                 0.03061849   FALSE  TRUE            FALSE
## PubDate.juliandate            1.39314146   FALSE FALSE            FALSE
## PubDate.zoo                  99.86221678   FALSE FALSE            FALSE
## UniqueID                    100.00000000   FALSE FALSE            FALSE
## .rnorm                      100.00000000   FALSE FALSE            FALSE
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
##          8   0.8944 0.57831   0.013981 0.06890         
##         16   0.9305 0.75939   0.004671 0.01720         
##         32   0.9303 0.75899   0.004681 0.01717         
##         55   0.9326 0.76882   0.004814 0.01705        *
## 
## The top 5 variables (out of 55):
##    WordCount.log1p, WordCount.root2, WordCount.nexp, NDSSName.my.fctrOpEd#Opnn#, PubDate.day.minutes.poly.1
## 
##  [1] "WordCount.log1p"                        
##  [2] "WordCount.root2"                        
##  [3] "WordCount.nexp"                         
##  [4] "NDSSName.my.fctrOpEd#Opnn#"             
##  [5] "PubDate.day.minutes.poly.1"             
##  [6] "PubDate.day.minutes.poly.4"             
##  [7] "PubDate.hour.fctr(15.3,23]"             
##  [8] "NDSSName.my.fctrScnc#Hlth#"             
##  [9] "PubDate.last4.log1p"                    
## [10] "PubDate.last2.log1p"                    
## [11] "NDSSName.my.fctrBsnss#Crsswrds/Gms#"    
## [12] "NDSSName.my.fctrStyls#U.S.#"            
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
## [26] "PubDate.minute.fctr(44.2,59.1]"         
## [27] "PubDate.day.minutes.poly.2"             
## [28] "PubDate.hour.fctr(7.67,15.3]"           
## [29] "PubDate.minute.fctr(14.8,29.5]"         
## [30] "PubDate.date.fctr(25,31]"               
## [31] "PubDate.last32.log1p"                   
## [32] "PubDate.second.fctr(44.2,59.1]"         
## [33] "PubDate.wkday.fctr3"                    
## [34] "NDSSName.my.fctrmyOthr"                 
## [35] "PubDate.date.fctr(19,25]"               
## [36] "NDSSName.my.fctr#Opnn#RmFrDbt"          
## [37] "NDSSName.my.fctrBsnss#Tchnlgy#"         
## [38] "PubDate.wkday.fctr4"                    
## [39] "PubDate.second.fctr(29.5,44.2]"         
## [40] "PubDate.date.fctr(13,19]"               
## [41] "NDSSName.my.fctrMtr#N.Y./Rgn#"          
## [42] "NDSSName.my.fctrTrvl#Trvl#"             
## [43] "NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss"
## [44] "NDSSName.my.fctr#Mltmd#"                
## [45] "PubDate.wkday.fctr2"                    
## [46] "NDSSName.my.fctrStyls##Fshn"            
## [47] "NDSSName.my.fctrFrgn#Wrld#"             
## [48] "PubDate.minute.fctr(29.5,44.2]"         
## [49] "NDSSName.my.fctrFrgn#Wrld#AsPcfc"       
## [50] "PubDate.wkday.fctr5"                    
## [51] "NDSSName.my.fctr#U.S.#Edctn"            
## [52] "NDSSName.my.fctrCltr#Arts#"             
## [53] "NDSSName.my.fctrBsnss#BsnssDy#Dlbk"     
## [54] "NDSSName.my.fctr##"                     
## [55] "NDSSName.my.fctrTStyl##"
```

![](NYTBlogs3_feat_PubDate_files/figure-html/select.features-3.png) 

```
## [1] "numeric data missing in : "
##   Popular Pplr.fctr 
##      1870      1870 
## [1] "numeric data w/ 0s in : "
##           WordCount             Popular     WordCount.log1p 
##                 109                5439                 109 
##     WordCount.root2      WordCount.nexp  PubDate.wkday.fctr 
##                 109                2044                 378 
##       PubDate.wkend       PubDate.hlday PubDate.day.minutes 
##                7787                8160                   5 
## [1] "numeric data w/ Infs in : "
## named integer(0)
## [1] "numeric data w/ NaNs in : "
## named integer(0)
## [1] "string data missing in : "
##       NewsDesk    SectionName SubsectionName       Headline        Snippet 
##           2408           2899           6176              0             13 
##       Abstract        PubDate    NDSSName.my           .lcn 
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
##              label step_major step_minor label_minor    bgn    end elapsed
## 9  select.features          5          0           0 62.061 83.208  21.147
## 10      fit.models          6          0           0 83.208     NA      NA
```

## Step `6.0: fit models`

```r
fit.models_0_chunk_df <- myadd_chunk(NULL, "fit.models_0_bgn", label.minor = "setup")
```

```
##              label step_major step_minor label_minor    bgn end elapsed
## 1 fit.models_0_bgn          1          0       setup 84.325  NA      NA
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
## 1 fit.models_0_bgn          1          0         setup 84.325 84.355
## 2 fit.models_0_MFO          1          1 myMFO_classfr 84.356     NA
##   elapsed
## 1   0.031
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
## 1 MFO###myMFO_classfr .rnorm               0                      0.296
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
## 2    fit.models_0_MFO          1          1    myMFO_classfr 84.356 87.427
## 3 fit.models_0_Random          1          2 myrandom_classfr 87.427     NA
##   elapsed
## 2   3.071
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
## 1                      0.298                 0.002       0.4990604
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
## 3 87.427 91.766   4.339
## 4 91.766     NA      NA
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
##                             (Intercept) 
##                             -4.57159198 
##                 NDSSName.my.fctr#Mltmd# 
##                             -1.22219085 
##           NDSSName.my.fctr#Opnn#RmFrDbt 
##                             -3.46072453 
##        NDSSName.my.fctr#Opnn#ThPblcEdtr 
##                              4.06871185 
##             NDSSName.my.fctr#U.S.#Edctn 
##                             -1.89443632 
##      NDSSName.my.fctrBsnss#BsnssDy#Dlbk 
##                             -0.22472818 
## NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss 
##                             -0.95537118 
##     NDSSName.my.fctrBsnss#Crsswrds/Gms# 
##                              4.55408513 
##          NDSSName.my.fctrBsnss#Tchnlgy# 
##                              0.77368538 
##              NDSSName.my.fctrCltr#Arts# 
##                             -0.09465691 
##              NDSSName.my.fctrFrgn#Wrld# 
##                             -1.45528874 
##        NDSSName.my.fctrFrgn#Wrld#AsPcfc 
##                             -1.60117505 
##           NDSSName.my.fctrMtr#N.Y./Rgn# 
##                              0.01563989 
##              NDSSName.my.fctrOpEd#Opnn# 
##                              4.51696382 
##              NDSSName.my.fctrScnc#Hlth# 
##                              3.51595317 
##             NDSSName.my.fctrStyls##Fshn 
##                             -1.85948925 
##             NDSSName.my.fctrStyls#U.S.# 
##                              3.27995325 
##                 NDSSName.my.fctrTStyl## 
##                             -1.54110404 
##              NDSSName.my.fctrTrvl#Trvl# 
##                             -1.41940605 
##                  NDSSName.my.fctrmyOthr 
##                             -1.90156922 
##                         WordCount.root2 
##                              0.08434378 
## [1] "max lambda < lambdaOpt:"
##                             (Intercept) 
##                             -4.60394059 
##                 NDSSName.my.fctr#Mltmd# 
##                             -1.25163328 
##           NDSSName.my.fctr#Opnn#RmFrDbt 
##                             -3.55521332 
##        NDSSName.my.fctr#Opnn#ThPblcEdtr 
##                              4.09217313 
##             NDSSName.my.fctr#U.S.#Edctn 
##                             -1.96172971 
##      NDSSName.my.fctrBsnss#BsnssDy#Dlbk 
##                             -0.22495986 
## NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss 
##                             -0.96836050 
##     NDSSName.my.fctrBsnss#Crsswrds/Gms# 
##                              4.58120497 
##          NDSSName.my.fctrBsnss#Tchnlgy# 
##                              0.78504703 
##              NDSSName.my.fctrCltr#Arts# 
##                             -0.09069661 
##              NDSSName.my.fctrFrgn#Wrld# 
##                             -1.51061232 
##        NDSSName.my.fctrFrgn#Wrld#AsPcfc 
##                             -1.63313235 
##           NDSSName.my.fctrMtr#N.Y./Rgn# 
##                              0.02466697 
##              NDSSName.my.fctrOpEd#Opnn# 
##                              4.54361134 
##              NDSSName.my.fctrScnc#Hlth# 
##                              3.53210055 
##             NDSSName.my.fctrStyls##Fshn 
##                             -1.92188290 
##             NDSSName.my.fctrStyls#U.S.# 
##                              3.29488750 
##                 NDSSName.my.fctrTStyl## 
##                             -1.57788931 
##              NDSSName.my.fctrTrvl#Trvl# 
##                             -1.47368131 
##                  NDSSName.my.fctrmyOthr 
##                             -1.97357582 
##                         WordCount.root2 
##                              0.08537319
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
##                           id                            feats
## 1 Max.cor.Y.rcv.1X1###glmnet WordCount.root2,NDSSName.my.fctr
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1               0                      1.024                 0.276
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
## [1] "fitting model: Max.cor.Y.rcv.3X1##rcv#glmnet"
## [1] "    indep_vars: WordCount.root2,NDSSName.my.fctr"
## Aggregating results
## Selecting tuning parameters
## Fitting alpha = 0.325, lambda = 0.0201 on full training set
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-12.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-13.png) 

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
##                             (Intercept) 
##                             -3.89350373 
##                 NDSSName.my.fctr#Mltmd# 
##                             -0.01916344 
##           NDSSName.my.fctr#Opnn#RmFrDbt 
##                             -1.18453357 
##        NDSSName.my.fctr#Opnn#ThPblcEdtr 
##                              3.21701058 
##             NDSSName.my.fctr#U.S.#Edctn 
##                             -0.47679040 
## NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss 
##                             -0.09891374 
##     NDSSName.my.fctrBsnss#Crsswrds/Gms# 
##                              3.87281404 
##          NDSSName.my.fctrBsnss#Tchnlgy# 
##                              0.40965256 
##              NDSSName.my.fctrFrgn#Wrld# 
##                             -0.05114617 
##        NDSSName.my.fctrFrgn#Wrld#AsPcfc 
##                             -0.47464340 
##              NDSSName.my.fctrOpEd#Opnn# 
##                              3.95214357 
##              NDSSName.my.fctrScnc#Hlth# 
##                              3.14232408 
##             NDSSName.my.fctrStyls##Fshn 
##                             -0.31867093 
##             NDSSName.my.fctrStyls#U.S.# 
##                              2.92610567 
##                 NDSSName.my.fctrTStyl## 
##                             -0.60538025 
##                         WordCount.root2 
##                              0.05783392 
## [1] "max lambda < lambdaOpt:"
##                             (Intercept) 
##                             -3.95644632 
##                 NDSSName.my.fctr#Mltmd# 
##                             -0.07182859 
##           NDSSName.my.fctr#Opnn#RmFrDbt 
##                             -1.30034382 
##        NDSSName.my.fctr#Opnn#ThPblcEdtr 
##                              3.29694236 
##             NDSSName.my.fctr#U.S.#Edctn 
##                             -0.53415905 
## NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss 
##                             -0.14259759 
##     NDSSName.my.fctrBsnss#Crsswrds/Gms# 
##                              3.94231812 
##          NDSSName.my.fctrBsnss#Tchnlgy# 
##                              0.45657914 
##              NDSSName.my.fctrFrgn#Wrld# 
##                             -0.10021084 
##        NDSSName.my.fctrFrgn#Wrld#AsPcfc 
##                             -0.53077048 
##              NDSSName.my.fctrOpEd#Opnn# 
##                              4.01324666 
##              NDSSName.my.fctrScnc#Hlth# 
##                              3.18936803 
##             NDSSName.my.fctrStyls##Fshn 
##                             -0.38069674 
##             NDSSName.my.fctrStyls#U.S.# 
##                              2.97176051 
##                 NDSSName.my.fctrTStyl## 
##                             -0.64837249 
##                         WordCount.root2 
##                              0.05978318
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-14.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-15.png) 

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

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-16.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-17.png) 

```
##          Prediction
## Reference    N    Y
##         N 1146  352
##         Y   67  163
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   7.575231e-01   3.107477e-01   7.365992e-01   7.775689e-01   8.668981e-01 
## AccuracyPValue  McnemarPValue 
##   1.000000e+00   9.066396e-44 
##                              id                            feats
## 1 Max.cor.Y.rcv.3X1##rcv#glmnet WordCount.root2,NDSSName.my.fctr
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1              25                      2.866                 0.286
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
## [1] "fitting model: Max.cor.Y.rcv.3X3##rcv#glmnet"
## [1] "    indep_vars: WordCount.root2,NDSSName.my.fctr"
## Aggregating results
## Selecting tuning parameters
## Fitting alpha = 0.325, lambda = 0.0201 on full training set
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-18.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-19.png) 

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
##                             (Intercept) 
##                             -3.89350373 
##                 NDSSName.my.fctr#Mltmd# 
##                             -0.01916344 
##           NDSSName.my.fctr#Opnn#RmFrDbt 
##                             -1.18453357 
##        NDSSName.my.fctr#Opnn#ThPblcEdtr 
##                              3.21701058 
##             NDSSName.my.fctr#U.S.#Edctn 
##                             -0.47679040 
## NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss 
##                             -0.09891374 
##     NDSSName.my.fctrBsnss#Crsswrds/Gms# 
##                              3.87281404 
##          NDSSName.my.fctrBsnss#Tchnlgy# 
##                              0.40965256 
##              NDSSName.my.fctrFrgn#Wrld# 
##                             -0.05114617 
##        NDSSName.my.fctrFrgn#Wrld#AsPcfc 
##                             -0.47464340 
##              NDSSName.my.fctrOpEd#Opnn# 
##                              3.95214357 
##              NDSSName.my.fctrScnc#Hlth# 
##                              3.14232408 
##             NDSSName.my.fctrStyls##Fshn 
##                             -0.31867093 
##             NDSSName.my.fctrStyls#U.S.# 
##                              2.92610567 
##                 NDSSName.my.fctrTStyl## 
##                             -0.60538025 
##                         WordCount.root2 
##                              0.05783392 
## [1] "max lambda < lambdaOpt:"
##                             (Intercept) 
##                             -3.95644632 
##                 NDSSName.my.fctr#Mltmd# 
##                             -0.07182859 
##           NDSSName.my.fctr#Opnn#RmFrDbt 
##                             -1.30034382 
##        NDSSName.my.fctr#Opnn#ThPblcEdtr 
##                              3.29694236 
##             NDSSName.my.fctr#U.S.#Edctn 
##                             -0.53415905 
## NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss 
##                             -0.14259759 
##     NDSSName.my.fctrBsnss#Crsswrds/Gms# 
##                              3.94231812 
##          NDSSName.my.fctrBsnss#Tchnlgy# 
##                              0.45657914 
##              NDSSName.my.fctrFrgn#Wrld# 
##                             -0.10021084 
##        NDSSName.my.fctrFrgn#Wrld#AsPcfc 
##                             -0.53077048 
##              NDSSName.my.fctrOpEd#Opnn# 
##                              4.01324666 
##              NDSSName.my.fctrScnc#Hlth# 
##                              3.18936803 
##             NDSSName.my.fctrStyls##Fshn 
##                             -0.38069674 
##             NDSSName.my.fctrStyls#U.S.# 
##                              2.97176051 
##                 NDSSName.my.fctrTStyl## 
##                             -0.64837249 
##                         WordCount.root2 
##                              0.05978318
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
##         N 1146  352
##         Y   67  163
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   7.575231e-01   3.107477e-01   7.365992e-01   7.775689e-01   8.668981e-01 
## AccuracyPValue  McnemarPValue 
##   1.000000e+00   9.066396e-44 
##                              id                            feats
## 1 Max.cor.Y.rcv.3X3##rcv#glmnet WordCount.root2,NDSSName.my.fctr
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1              25                      4.855                 0.273
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
## [1] "fitting model: Max.cor.Y.rcv.3X5##rcv#glmnet"
## [1] "    indep_vars: WordCount.root2,NDSSName.my.fctr"
## Aggregating results
## Selecting tuning parameters
## Fitting alpha = 0.325, lambda = 0.0201 on full training set
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-24.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-25.png) 

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
##                             (Intercept) 
##                             -3.89350373 
##                 NDSSName.my.fctr#Mltmd# 
##                             -0.01916344 
##           NDSSName.my.fctr#Opnn#RmFrDbt 
##                             -1.18453357 
##        NDSSName.my.fctr#Opnn#ThPblcEdtr 
##                              3.21701058 
##             NDSSName.my.fctr#U.S.#Edctn 
##                             -0.47679040 
## NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss 
##                             -0.09891374 
##     NDSSName.my.fctrBsnss#Crsswrds/Gms# 
##                              3.87281404 
##          NDSSName.my.fctrBsnss#Tchnlgy# 
##                              0.40965256 
##              NDSSName.my.fctrFrgn#Wrld# 
##                             -0.05114617 
##        NDSSName.my.fctrFrgn#Wrld#AsPcfc 
##                             -0.47464340 
##              NDSSName.my.fctrOpEd#Opnn# 
##                              3.95214357 
##              NDSSName.my.fctrScnc#Hlth# 
##                              3.14232408 
##             NDSSName.my.fctrStyls##Fshn 
##                             -0.31867093 
##             NDSSName.my.fctrStyls#U.S.# 
##                              2.92610567 
##                 NDSSName.my.fctrTStyl## 
##                             -0.60538025 
##                         WordCount.root2 
##                              0.05783392 
## [1] "max lambda < lambdaOpt:"
##                             (Intercept) 
##                             -3.95644632 
##                 NDSSName.my.fctr#Mltmd# 
##                             -0.07182859 
##           NDSSName.my.fctr#Opnn#RmFrDbt 
##                             -1.30034382 
##        NDSSName.my.fctr#Opnn#ThPblcEdtr 
##                              3.29694236 
##             NDSSName.my.fctr#U.S.#Edctn 
##                             -0.53415905 
## NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss 
##                             -0.14259759 
##     NDSSName.my.fctrBsnss#Crsswrds/Gms# 
##                              3.94231812 
##          NDSSName.my.fctrBsnss#Tchnlgy# 
##                              0.45657914 
##              NDSSName.my.fctrFrgn#Wrld# 
##                             -0.10021084 
##        NDSSName.my.fctrFrgn#Wrld#AsPcfc 
##                             -0.53077048 
##              NDSSName.my.fctrOpEd#Opnn# 
##                              4.01324666 
##              NDSSName.my.fctrScnc#Hlth# 
##                              3.18936803 
##             NDSSName.my.fctrStyls##Fshn 
##                             -0.38069674 
##             NDSSName.my.fctrStyls#U.S.# 
##                              2.97176051 
##                 NDSSName.my.fctrTStyl## 
##                             -0.64837249 
##                         WordCount.root2 
##                              0.05978318
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-26.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-27.png) 

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

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-28.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-29.png) 

```
##          Prediction
## Reference    N    Y
##         N 1146  352
##         Y   67  163
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   7.575231e-01   3.107477e-01   7.365992e-01   7.775689e-01   8.668981e-01 
## AccuracyPValue  McnemarPValue 
##   1.000000e+00   9.066396e-44 
##                              id                            feats
## 1 Max.cor.Y.rcv.3X5##rcv#glmnet WordCount.root2,NDSSName.my.fctr
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1              25                      8.311                 0.277
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
## [1] "fitting model: Max.cor.Y.rcv.5X1##rcv#glmnet"
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

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-30.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-31.png) 

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
##                             (Intercept) 
##                             -3.81141260 
##                 NDSSName.my.fctr#Mltmd# 
##                             -0.68105584 
##           NDSSName.my.fctr#Opnn#RmFrDbt 
##                             -1.92624537 
##        NDSSName.my.fctr#Opnn#ThPblcEdtr 
##                              3.40699589 
##             NDSSName.my.fctr#U.S.#Edctn 
##                             -0.98291999 
##      NDSSName.my.fctrBsnss#BsnssDy#Dlbk 
##                             -0.22577146 
## NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss 
##                             -0.64343834 
##     NDSSName.my.fctrBsnss#Crsswrds/Gms# 
##                              3.82797332 
##          NDSSName.my.fctrBsnss#Tchnlgy# 
##                              0.45317927 
##              NDSSName.my.fctrCltr#Arts# 
##                             -0.17187706 
##              NDSSName.my.fctrFrgn#Wrld# 
##                             -0.72035867 
##        NDSSName.my.fctrFrgn#Wrld#AsPcfc 
##                             -0.99018968 
##              NDSSName.my.fctrOpEd#Opnn# 
##                              3.81891156 
##              NDSSName.my.fctrScnc#Hlth# 
##                              3.05516080 
##             NDSSName.my.fctrStyls##Fshn 
##                             -0.97651721 
##             NDSSName.my.fctrStyls#U.S.# 
##                              2.84779285 
##                 NDSSName.my.fctrTStyl## 
##                             -0.94109645 
##              NDSSName.my.fctrTrvl#Trvl# 
##                             -0.68827560 
##                  NDSSName.my.fctrmyOthr 
##                             -0.84423735 
##                         WordCount.root2 
##                              0.06115867 
## [1] "max lambda < lambdaOpt:"
##                             (Intercept) 
##                             -3.87108412 
##                 NDSSName.my.fctr#Mltmd# 
##                             -0.71588942 
##           NDSSName.my.fctr#Opnn#RmFrDbt 
##                             -2.02010163 
##        NDSSName.my.fctr#Opnn#ThPblcEdtr 
##                              3.46715540 
##             NDSSName.my.fctr#U.S.#Edctn 
##                             -1.02957582 
##      NDSSName.my.fctrBsnss#BsnssDy#Dlbk 
##                             -0.22558850 
## NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss 
##                             -0.66798026 
##     NDSSName.my.fctrBsnss#Crsswrds/Gms# 
##                              3.89132347 
##          NDSSName.my.fctrBsnss#Tchnlgy# 
##                              0.48212450 
##              NDSSName.my.fctrCltr#Arts# 
##                             -0.16733777 
##              NDSSName.my.fctrFrgn#Wrld# 
##                             -0.75793881 
##        NDSSName.my.fctrFrgn#Wrld#AsPcfc 
##                             -1.03076807 
##              NDSSName.my.fctrOpEd#Opnn# 
##                              3.87908175 
##              NDSSName.my.fctrScnc#Hlth# 
##                              3.09788786 
##             NDSSName.my.fctrStyls##Fshn 
##                             -1.02481879 
##             NDSSName.my.fctrStyls#U.S.# 
##                              2.88826078 
##                 NDSSName.my.fctrTStyl## 
##                             -0.97585470 
##              NDSSName.my.fctrTrvl#Trvl# 
##                             -0.72668427 
##                  NDSSName.my.fctrmyOthr 
##                             -0.90347045 
##                         WordCount.root2 
##                              0.06289698
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-32.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-33.png) 

```
##          Prediction
## Reference    N    Y
##         N 3800  141
##         Y  179  684
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   9.333888e-01   7.700473e-01   9.259666e-01   9.402789e-01   8.203580e-01 
## AccuracyPValue  McnemarPValue 
##  1.097051e-115   3.860591e-02
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-34.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-35.png) 

```
##          Prediction
## Reference    N    Y
##         N 1137  361
##         Y   53  177
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   7.604167e-01   3.373693e-01   7.395703e-01   7.803749e-01   8.668981e-01 
## AccuracyPValue  McnemarPValue 
##   1.000000e+00   1.935747e-51 
##                              id                            feats
## 1 Max.cor.Y.rcv.5X1##rcv#glmnet WordCount.root2,NDSSName.my.fctr
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1              25                      3.432                 0.271
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
## [1] "fitting model: Max.cor.Y.rcv.5X3##rcv#glmnet"
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

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-36.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-37.png) 

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
##                             (Intercept) 
##                             -3.81141260 
##                 NDSSName.my.fctr#Mltmd# 
##                             -0.68105584 
##           NDSSName.my.fctr#Opnn#RmFrDbt 
##                             -1.92624537 
##        NDSSName.my.fctr#Opnn#ThPblcEdtr 
##                              3.40699589 
##             NDSSName.my.fctr#U.S.#Edctn 
##                             -0.98291999 
##      NDSSName.my.fctrBsnss#BsnssDy#Dlbk 
##                             -0.22577146 
## NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss 
##                             -0.64343834 
##     NDSSName.my.fctrBsnss#Crsswrds/Gms# 
##                              3.82797332 
##          NDSSName.my.fctrBsnss#Tchnlgy# 
##                              0.45317927 
##              NDSSName.my.fctrCltr#Arts# 
##                             -0.17187706 
##              NDSSName.my.fctrFrgn#Wrld# 
##                             -0.72035867 
##        NDSSName.my.fctrFrgn#Wrld#AsPcfc 
##                             -0.99018968 
##              NDSSName.my.fctrOpEd#Opnn# 
##                              3.81891156 
##              NDSSName.my.fctrScnc#Hlth# 
##                              3.05516080 
##             NDSSName.my.fctrStyls##Fshn 
##                             -0.97651721 
##             NDSSName.my.fctrStyls#U.S.# 
##                              2.84779285 
##                 NDSSName.my.fctrTStyl## 
##                             -0.94109645 
##              NDSSName.my.fctrTrvl#Trvl# 
##                             -0.68827560 
##                  NDSSName.my.fctrmyOthr 
##                             -0.84423735 
##                         WordCount.root2 
##                              0.06115867 
## [1] "max lambda < lambdaOpt:"
##                             (Intercept) 
##                             -3.87108412 
##                 NDSSName.my.fctr#Mltmd# 
##                             -0.71588942 
##           NDSSName.my.fctr#Opnn#RmFrDbt 
##                             -2.02010163 
##        NDSSName.my.fctr#Opnn#ThPblcEdtr 
##                              3.46715540 
##             NDSSName.my.fctr#U.S.#Edctn 
##                             -1.02957582 
##      NDSSName.my.fctrBsnss#BsnssDy#Dlbk 
##                             -0.22558850 
## NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss 
##                             -0.66798026 
##     NDSSName.my.fctrBsnss#Crsswrds/Gms# 
##                              3.89132347 
##          NDSSName.my.fctrBsnss#Tchnlgy# 
##                              0.48212450 
##              NDSSName.my.fctrCltr#Arts# 
##                             -0.16733777 
##              NDSSName.my.fctrFrgn#Wrld# 
##                             -0.75793881 
##        NDSSName.my.fctrFrgn#Wrld#AsPcfc 
##                             -1.03076807 
##              NDSSName.my.fctrOpEd#Opnn# 
##                              3.87908175 
##              NDSSName.my.fctrScnc#Hlth# 
##                              3.09788786 
##             NDSSName.my.fctrStyls##Fshn 
##                             -1.02481879 
##             NDSSName.my.fctrStyls#U.S.# 
##                              2.88826078 
##                 NDSSName.my.fctrTStyl## 
##                             -0.97585470 
##              NDSSName.my.fctrTrvl#Trvl# 
##                             -0.72668427 
##                  NDSSName.my.fctrmyOthr 
##                             -0.90347045 
##                         WordCount.root2 
##                              0.06289698
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-38.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-39.png) 

```
##          Prediction
## Reference    N    Y
##         N 3800  141
##         Y  179  684
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   9.333888e-01   7.700473e-01   9.259666e-01   9.402789e-01   8.203580e-01 
## AccuracyPValue  McnemarPValue 
##  1.097051e-115   3.860591e-02
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-40.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-41.png) 

```
##          Prediction
## Reference    N    Y
##         N 1137  361
##         Y   53  177
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   7.604167e-01   3.373693e-01   7.395703e-01   7.803749e-01   8.668981e-01 
## AccuracyPValue  McnemarPValue 
##   1.000000e+00   1.935747e-51 
##                              id                            feats
## 1 Max.cor.Y.rcv.5X3##rcv#glmnet WordCount.root2,NDSSName.my.fctr
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1              25                      6.452                 0.277
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
## [1] "fitting model: Max.cor.Y.rcv.5X5##rcv#glmnet"
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

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-42.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-43.png) 

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
##                             (Intercept) 
##                             -3.81141260 
##                 NDSSName.my.fctr#Mltmd# 
##                             -0.68105584 
##           NDSSName.my.fctr#Opnn#RmFrDbt 
##                             -1.92624537 
##        NDSSName.my.fctr#Opnn#ThPblcEdtr 
##                              3.40699589 
##             NDSSName.my.fctr#U.S.#Edctn 
##                             -0.98291999 
##      NDSSName.my.fctrBsnss#BsnssDy#Dlbk 
##                             -0.22577146 
## NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss 
##                             -0.64343834 
##     NDSSName.my.fctrBsnss#Crsswrds/Gms# 
##                              3.82797332 
##          NDSSName.my.fctrBsnss#Tchnlgy# 
##                              0.45317927 
##              NDSSName.my.fctrCltr#Arts# 
##                             -0.17187706 
##              NDSSName.my.fctrFrgn#Wrld# 
##                             -0.72035867 
##        NDSSName.my.fctrFrgn#Wrld#AsPcfc 
##                             -0.99018968 
##              NDSSName.my.fctrOpEd#Opnn# 
##                              3.81891156 
##              NDSSName.my.fctrScnc#Hlth# 
##                              3.05516080 
##             NDSSName.my.fctrStyls##Fshn 
##                             -0.97651721 
##             NDSSName.my.fctrStyls#U.S.# 
##                              2.84779285 
##                 NDSSName.my.fctrTStyl## 
##                             -0.94109645 
##              NDSSName.my.fctrTrvl#Trvl# 
##                             -0.68827560 
##                  NDSSName.my.fctrmyOthr 
##                             -0.84423735 
##                         WordCount.root2 
##                              0.06115867 
## [1] "max lambda < lambdaOpt:"
##                             (Intercept) 
##                             -3.87108412 
##                 NDSSName.my.fctr#Mltmd# 
##                             -0.71588942 
##           NDSSName.my.fctr#Opnn#RmFrDbt 
##                             -2.02010163 
##        NDSSName.my.fctr#Opnn#ThPblcEdtr 
##                              3.46715540 
##             NDSSName.my.fctr#U.S.#Edctn 
##                             -1.02957582 
##      NDSSName.my.fctrBsnss#BsnssDy#Dlbk 
##                             -0.22558850 
## NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss 
##                             -0.66798026 
##     NDSSName.my.fctrBsnss#Crsswrds/Gms# 
##                              3.89132347 
##          NDSSName.my.fctrBsnss#Tchnlgy# 
##                              0.48212450 
##              NDSSName.my.fctrCltr#Arts# 
##                             -0.16733777 
##              NDSSName.my.fctrFrgn#Wrld# 
##                             -0.75793881 
##        NDSSName.my.fctrFrgn#Wrld#AsPcfc 
##                             -1.03076807 
##              NDSSName.my.fctrOpEd#Opnn# 
##                              3.87908175 
##              NDSSName.my.fctrScnc#Hlth# 
##                              3.09788786 
##             NDSSName.my.fctrStyls##Fshn 
##                             -1.02481879 
##             NDSSName.my.fctrStyls#U.S.# 
##                              2.88826078 
##                 NDSSName.my.fctrTStyl## 
##                             -0.97585470 
##              NDSSName.my.fctrTrvl#Trvl# 
##                             -0.72668427 
##                  NDSSName.my.fctrmyOthr 
##                             -0.90347045 
##                         WordCount.root2 
##                              0.06289698
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-44.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-45.png) 

```
##          Prediction
## Reference    N    Y
##         N 3800  141
##         Y  179  684
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   9.333888e-01   7.700473e-01   9.259666e-01   9.402789e-01   8.203580e-01 
## AccuracyPValue  McnemarPValue 
##  1.097051e-115   3.860591e-02
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-46.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-47.png) 

```
##          Prediction
## Reference    N    Y
##         N 1137  361
##         Y   53  177
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   7.604167e-01   3.373693e-01   7.395703e-01   7.803749e-01   8.668981e-01 
## AccuracyPValue  McnemarPValue 
##   1.000000e+00   1.935747e-51 
##                              id                            feats
## 1 Max.cor.Y.rcv.5X5##rcv#glmnet WordCount.root2,NDSSName.my.fctr
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1              25                       9.15                 0.269
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

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-48.png) 

```r
fit.models_0_chunk_df <- myadd_chunk(fit.models_0_chunk_df, 
                    paste0("fit.models_0_", "Max.cor.Y[rcv.1X1.cp.0|]"), major.inc = FALSE,
                                    label.minor = "rpart")
```

```
##                                   label step_major step_minor label_minor
## 4        fit.models_0_Max.cor.Y.rcv.*X*          1          3      glmnet
## 5 fit.models_0_Max.cor.Y[rcv.1X1.cp.0|]          1          4       rpart
##       bgn     end elapsed
## 4  91.766 169.327  77.561
## 5 169.327      NA      NA
```

```r
ret_lst <- myfit_mdl(mdl_specs_lst=myinit_mdl_specs_lst(mdl_specs_lst=list(
    id.prefix="Max.cor.Y.rcv.1X1.cp.0", type=glb_model_type, trainControl.method="none",
    train.method="rpart",
    tune.df=data.frame(method="rpart", parameter="cp", min=0.0, max=0.0, by=0.1))),
                    indep_vars=max_cor_y_x_vars, rsp_var=glb_rsp_var, 
                    fit_df=glbObsFit, OOB_df=glbObsOOB)
```

```
## [1] "fitting model: Max.cor.Y.rcv.1X1.cp.0###rpart"
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

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-49.png) 

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
##          NDSSName.my.fctrOpEd#Opnn# NDSSName.my.fctrBsnss#Crsswrds/Gms# 
##                                  48                                  14 
##          NDSSName.my.fctrScnc#Hlth#         NDSSName.my.fctrStyls#U.S.# 
##                                  14                                  11 
##                     WordCount.root2    NDSSName.my.fctr#Opnn#ThPblcEdtr 
##                                   9                                   2 
##  NDSSName.my.fctrBsnss#BsnssDy#Dlbk 
##                                   1 
## 
## Node number 1: 4804 observations,    complexity param=0.3696408
##   predicted class=N  expected loss=0.179642  P(node) =1
##     class counts:  3941   863
##    probabilities: 0.820 0.180 
##   left son=2 (4367 obs) right son=3 (437 obs)
##   Primary splits:
##       NDSSName.my.fctrOpEd#Opnn#          < 0.5      to the left,  improve=451.59770, (0 missing)
##       NDSSName.my.fctrBsnss#Crsswrds/Gms# < 0.5      to the left,  improve=112.88510, (0 missing)
##       WordCount.root2                     < 25.75849 to the left,  improve=111.17610, (0 missing)
##       NDSSName.my.fctrScnc#Hlth#          < 0.5      to the left,  improve= 99.35206, (0 missing)
##       NDSSName.my.fctrStyls#U.S.#         < 0.5      to the left,  improve= 68.73272, (0 missing)
## 
## Node number 2: 4367 observations,    complexity param=0.09849363
##   predicted class=N  expected loss=0.1110602  P(node) =0.9090341
##     class counts:  3882   485
##    probabilities: 0.889 0.111 
##   left son=4 (4262 obs) right son=5 (105 obs)
##   Primary splits:
##       NDSSName.my.fctrBsnss#Crsswrds/Gms# < 0.5      to the left,  improve=135.55130, (0 missing)
##       NDSSName.my.fctrScnc#Hlth#          < 0.5      to the left,  improve=125.07920, (0 missing)
##       WordCount.root2                     < 25.75849 to the left,  improve= 94.70710, (0 missing)
##       NDSSName.my.fctrStyls#U.S.#         < 0.5      to the left,  improve= 88.56821, (0 missing)
##       NDSSName.my.fctr#Opnn#ThPblcEdtr    < 0.5      to the left,  improve= 18.74400, (0 missing)
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
##       NDSSName.my.fctrScnc#Hlth#       < 0.5      to the left,  improve=132.96710, (0 missing)
##       NDSSName.my.fctrStyls#U.S.#      < 0.5      to the left,  improve= 94.69099, (0 missing)
##       WordCount.root2                  < 26.49528 to the left,  improve= 84.07487, (0 missing)
##       NDSSName.my.fctr#Opnn#ThPblcEdtr < 0.5      to the left,  improve= 19.71762, (0 missing)
##       NDSSName.my.fctrTStyl##          < 0.5      to the right, improve= 10.17000, (0 missing)
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
##       NDSSName.my.fctrStyls#U.S.#      < 0.5      to the left,  improve=102.410700, (0 missing)
##       WordCount.root2                  < 25.01    to the left,  improve= 47.352210, (0 missing)
##       NDSSName.my.fctr#Opnn#ThPblcEdtr < 0.5      to the left,  improve= 20.930810, (0 missing)
##       NDSSName.my.fctrTStyl##          < 0.5      to the right, improve=  5.249425, (0 missing)
##       NDSSName.my.fctrBsnss#Tchnlgy#   < 0.5      to the left,  improve=  2.395935, (0 missing)
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
##       WordCount.root2                  < 25.01    to the left,  improve=29.253580, (0 missing)
##       NDSSName.my.fctr#Opnn#ThPblcEdtr < 0.5      to the left,  improve=21.978920, (0 missing)
##       NDSSName.my.fctrBsnss#Tchnlgy#   < 0.5      to the left,  improve= 3.887348, (0 missing)
##       NDSSName.my.fctrTStyl##          < 0.5      to the right, improve= 2.348653, (0 missing)
##       NDSSName.my.fctr#U.S.#Edctn      < 0.5      to the right, improve= 1.187739, (0 missing)
##   Surrogate splits:
##       NDSSName.my.fctr#Opnn#RmFrDbt           < 0.5      to the left,  agree=0.758, adj=0.042, (0 split)
##       NDSSName.my.fctrFrgn#Wrld#AsPcfc        < 0.5      to the left,  agree=0.752, adj=0.016, (0 split)
##       NDSSName.my.fctr#Opnn#ThPblcEdtr        < 0.5      to the left,  agree=0.750, adj=0.008, (0 split)
##       NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss < 0.5      to the left,  agree=0.748, adj=0.002, (0 split)
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
##       NDSSName.my.fctr#Opnn#ThPblcEdtr < 0.5      to the left,  improve=14.193880, (0 missing)
##       NDSSName.my.fctrCltr#Arts#       < 0.5      to the left,  improve= 3.669601, (0 missing)
##       NDSSName.my.fctrFrgn#Wrld#AsPcfc < 0.5      to the right, improve= 3.556158, (0 missing)
##       WordCount.root2                  < 34.19795 to the left,  improve= 2.582851, (0 missing)
##       NDSSName.my.fctr#Opnn#RmFrDbt    < 0.5      to the right, improve= 2.031748, (0 missing)
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
##       NDSSName.my.fctrCltr#Arts#       < 0.5      to the left,  improve=4.094729, (0 missing)
##       NDSSName.my.fctrFrgn#Wrld#AsPcfc < 0.5      to the right, improve=3.106316, (0 missing)
##       WordCount.root2                  < 29.5127  to the left,  improve=2.722793, (0 missing)
##       NDSSName.my.fctrBsnss#Tchnlgy#   < 0.5      to the left,  improve=1.962300, (0 missing)
##       NDSSName.my.fctr#Opnn#RmFrDbt    < 0.5      to the right, improve=1.793603, (0 missing)
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
##       WordCount.root2                    < 33.97057 to the left,  improve=2.913816, (0 missing)
##       NDSSName.my.fctrFrgn#Wrld#AsPcfc   < 0.5      to the right, improve=2.586923, (0 missing)
##       NDSSName.my.fctrBsnss#Tchnlgy#     < 0.5      to the left,  improve=2.402029, (0 missing)
##       NDSSName.my.fctr#Opnn#RmFrDbt      < 0.5      to the right, improve=1.513920, (0 missing)
##       NDSSName.my.fctrBsnss#BsnssDy#Dlbk < 0.5      to the left,  improve=1.276783, (0 missing)
##   Surrogate splits:
##       NDSSName.my.fctr#Opnn#RmFrDbt < 0.5      to the left,  agree=0.719, adj=0.139, (0 split)
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
##       NDSSName.my.fctrBsnss#Tchnlgy#   < 0.5      to the left,  improve=2.8404170, (0 missing)
##       NDSSName.my.fctrFrgn#Wrld#AsPcfc < 0.5      to the right, improve=1.0796950, (0 missing)
##       NDSSName.my.fctrTStyl##          < 0.5      to the right, improve=1.0670160, (0 missing)
##       WordCount.root2                  < 29.5127  to the left,  improve=0.8966879, (0 missing)
##       NDSSName.my.fctr#Mltmd#          < 0.5      to the right, improve=0.4399337, (0 missing)
## 
## Node number 265: 303 observations,    complexity param=0.0005793743
##   predicted class=N  expected loss=0.1881188  P(node) =0.06307244
##     class counts:   246    57
##    probabilities: 0.812 0.188 
##   left son=530 (222 obs) right son=531 (81 obs)
##   Primary splits:
##       NDSSName.my.fctrBsnss#BsnssDy#Dlbk < 0.5      to the left,  improve=5.4890570, (0 missing)
##       WordCount.root2                    < 38.17067 to the right, improve=5.0156320, (0 missing)
##       NDSSName.my.fctr#Opnn#RmFrDbt      < 0.5      to the right, improve=3.4510070, (0 missing)
##       NDSSName.my.fctrFrgn#Wrld#AsPcfc   < 0.5      to the right, improve=1.5155860, (0 missing)
##       NDSSName.my.fctr#U.S.#Edctn        < 0.5      to the right, improve=0.8078801, (0 missing)
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
##       WordCount.root2                  < 29.33428 to the left,  improve=1.5853030, (0 missing)
##       NDSSName.my.fctrTStyl##          < 0.5      to the right, improve=0.7645570, (0 missing)
##       NDSSName.my.fctrFrgn#Wrld#AsPcfc < 0.5      to the right, improve=0.7250433, (0 missing)
##       NDSSName.my.fctrStyls##Fshn      < 0.5      to the right, improve=0.3000638, (0 missing)
##       NDSSName.my.fctr#Mltmd#          < 0.5      to the right, improve=0.2729836, (0 missing)
##   Surrogate splits:
##       NDSSName.my.fctrMtr#N.Y./Rgn#           < 0.5      to the left,  agree=0.560, adj=0.118, (0 split)
##       NDSSName.my.fctrFrgn#Wrld#AsPcfc        < 0.5      to the right, agree=0.533, adj=0.064, (0 split)
##       NDSSName.my.fctr#Mltmd#                 < 0.5      to the right, agree=0.524, adj=0.046, (0 split)
##       NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss < 0.5      to the left,  agree=0.515, adj=0.029, (0 split)
##       NDSSName.my.fctrStyls##Fshn             < 0.5      to the right, agree=0.512, adj=0.021, (0 split)
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
##       WordCount.root2                  < 32.57299 to the right, improve=0.8968765, (0 missing)
##       NDSSName.my.fctrTStyl##          < 0.5      to the right, improve=0.7830739, (0 missing)
##       NDSSName.my.fctrMtr#N.Y./Rgn#    < 0.5      to the right, improve=0.3683673, (0 missing)
##       NDSSName.my.fctr#Mltmd#          < 0.5      to the right, improve=0.3578067, (0 missing)
##       NDSSName.my.fctrFrgn#Wrld#AsPcfc < 0.5      to the right, improve=0.3021494, (0 missing)
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
##       NDSSName.my.fctrTStyl##          < 0.5      to the right, improve=0.5601729, (0 missing)
##       NDSSName.my.fctr#Mltmd#          < 0.5      to the right, improve=0.5108985, (0 missing)
##       WordCount.root2                  < 30.09153 to the right, improve=0.4980706, (0 missing)
##       NDSSName.my.fctrMtr#N.Y./Rgn#    < 0.5      to the right, improve=0.4241343, (0 missing)
##       NDSSName.my.fctrFrgn#Wrld#AsPcfc < 0.5      to the right, improve=0.3390226, (0 missing)
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
##       NDSSName.my.fctr#Mltmd#                 < 0.5      to the right, improve=0.5769882, (0 missing)
##       NDSSName.my.fctrMtr#N.Y./Rgn#           < 0.5      to the right, improve=0.5314217, (0 missing)
##       WordCount.root2                         < 30.09153 to the right, improve=0.4682049, (0 missing)
##       NDSSName.my.fctrFrgn#Wrld#AsPcfc        < 0.5      to the right, improve=0.4106319, (0 missing)
##       NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss < 0.5      to the right, improve=0.1814254, (0 missing)
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
##       NDSSName.my.fctrMtr#N.Y./Rgn#           < 0.5      to the right, improve=0.6559045, (0 missing)
##       NDSSName.my.fctrFrgn#Wrld#AsPcfc        < 0.5      to the right, improve=0.4920635, (0 missing)
##       WordCount.root2                         < 30.09153 to the right, improve=0.3890196, (0 missing)
##       NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss < 0.5      to the right, improve=0.2415584, (0 missing)
##       NDSSName.my.fctrBsnss#BsnssDy#Dlbk      < 0.5      to the left,  improve=0.0126479, (0 missing)
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
##       NDSSName.my.fctrFrgn#Wrld#AsPcfc        < 0.5      to the right, improve=0.67831090, (0 missing)
##       WordCount.root2                         < 32.38827 to the left,  improve=0.61044970, (0 missing)
##       NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss < 0.5      to the right, improve=0.38816480, (0 missing)
##       NDSSName.my.fctrBsnss#BsnssDy#Dlbk      < 0.5      to the right, improve=0.01539613, (0 missing)
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
##       WordCount.root2                         < 30.09153 to the right, improve=0.9266317, (0 missing)
##       NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss < 0.5      to the right, improve=0.5580040, (0 missing)
##       NDSSName.my.fctrBsnss#BsnssDy#Dlbk      < 0.5      to the right, improve=0.1306354, (0 missing)
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
##       WordCount.root2                    < 29.92488 to the left,  improve=3.00231700, (0 missing)
##       NDSSName.my.fctrBsnss#BsnssDy#Dlbk < 0.5      to the left,  improve=0.01303089, (0 missing)
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
##        2) NDSSName.my.fctrOpEd#Opnn#< 0.5 4367 485 N (0.88893978 0.11106022)  
##          4) NDSSName.my.fctrBsnss#Crsswrds/Gms#< 0.5 4262 390 N (0.90849366 0.09150634)  
##            8) NDSSName.my.fctrScnc#Hlth#< 0.5 4114 279 N (0.93218279 0.06781721)  
##             16) NDSSName.my.fctrStyls#U.S.#< 0.5 3987 191 N (0.95209431 0.04790569)  
##               32) WordCount.root2< 25.01 2982  38 N (0.98725687 0.01274313) *
##               33) WordCount.root2>=25.01 1005 153 N (0.84776119 0.15223881)  
##                 66) NDSSName.my.fctr#Opnn#ThPblcEdtr< 0.5 993 142 N (0.85699899 0.14300101)  
##                  132) NDSSName.my.fctrCltr#Arts#< 0.5 930 122 N (0.86881720 0.13118280)  
##                    264) WordCount.root2< 33.97057 627  65 N (0.89633174 0.10366826)  
##                      528) NDSSName.my.fctrBsnss#Tchnlgy#< 0.5 561  49 N (0.91265597 0.08734403)  
##                       1056) WordCount.root2< 29.33428 281  14 N (0.95017794 0.04982206) *
##                       1057) WordCount.root2>=29.33428 280  35 N (0.87500000 0.12500000)  
##                         2114) WordCount.root2>=32.57299 71   4 N (0.94366197 0.05633803) *
##                         2115) WordCount.root2< 32.57299 209  31 N (0.85167464 0.14832536)  
##                           4230) NDSSName.my.fctrTStyl##>=0.5 12   0 N (1.00000000 0.00000000) *
##                           4231) NDSSName.my.fctrTStyl##< 0.5 197  31 N (0.84263959 0.15736041)  
##                             8462) NDSSName.my.fctr#Mltmd#>=0.5 11   0 N (1.00000000 0.00000000) *
##                             8463) NDSSName.my.fctr#Mltmd#< 0.5 186  31 N (0.83333333 0.16666667)  
##                              16926) NDSSName.my.fctrMtr#N.Y./Rgn#>=0.5 29   2 N (0.93103448 0.06896552) *
##                              16927) NDSSName.my.fctrMtr#N.Y./Rgn#< 0.5 157  29 N (0.81528662 0.18471338)  
##                                33854) NDSSName.my.fctrFrgn#Wrld#AsPcfc>=0.5 18   1 N (0.94444444 0.05555556) *
##                                33855) NDSSName.my.fctrFrgn#Wrld#AsPcfc< 0.5 139  28 N (0.79856115 0.20143885)  
##                                  67710) WordCount.root2>=30.09153 102  17 N (0.83333333 0.16666667) *
##                                  67711) WordCount.root2< 30.09153 37  11 N (0.70270270 0.29729730)  
##                                   135422) WordCount.root2< 29.92488 30   6 N (0.80000000 0.20000000) *
##                                   135423) WordCount.root2>=29.92488 7   2 Y (0.28571429 0.71428571) *
##                      529) NDSSName.my.fctrBsnss#Tchnlgy#>=0.5 66  16 N (0.75757576 0.24242424)  
##                       1058) WordCount.root2< 27.86575 38   7 N (0.81578947 0.18421053) *
##                       1059) WordCount.root2>=27.86575 28   9 N (0.67857143 0.32142857)  
##                         2118) WordCount.root2>=28.6269 19   4 N (0.78947368 0.21052632) *
##                         2119) WordCount.root2< 28.6269 9   4 Y (0.44444444 0.55555556) *
##                    265) WordCount.root2>=33.97057 303  57 N (0.81188119 0.18811881)  
##                      530) NDSSName.my.fctrBsnss#BsnssDy#Dlbk< 0.5 222  29 N (0.86936937 0.13063063) *
##                      531) NDSSName.my.fctrBsnss#BsnssDy#Dlbk>=0.5 81  28 N (0.65432099 0.34567901)  
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
##                  133) NDSSName.my.fctrCltr#Arts#>=0.5 63  20 N (0.68253968 0.31746032)  
##                    266) WordCount.root2< 26.99984 14   3 N (0.78571429 0.21428571) *
##                    267) WordCount.root2>=26.99984 49  17 N (0.65306122 0.34693878)  
##                      534) WordCount.root2>=41.56249 10   2 N (0.80000000 0.20000000) *
##                      535) WordCount.root2< 41.56249 39  15 N (0.61538462 0.38461538)  
##                       1070) WordCount.root2< 34.23387 32  11 N (0.65625000 0.34375000) *
##                       1071) WordCount.root2>=34.23387 7   3 Y (0.42857143 0.57142857) *
##                 67) NDSSName.my.fctr#Opnn#ThPblcEdtr>=0.5 12   1 Y (0.08333333 0.91666667) *
##             17) NDSSName.my.fctrStyls#U.S.#>=0.5 127  39 Y (0.30708661 0.69291339)  
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
##            9) NDSSName.my.fctrScnc#Hlth#>=0.5 148  37 Y (0.25000000 0.75000000)  
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
##          5) NDSSName.my.fctrBsnss#Crsswrds/Gms#>=0.5 105  10 Y (0.09523810 0.90476190)  
##           10) WordCount.root2< 18.9043 12   5 N (0.58333333 0.41666667) *
##           11) WordCount.root2>=18.9043 93   3 Y (0.03225806 0.96774194) *
##        3) NDSSName.my.fctrOpEd#Opnn#>=0.5 437  59 Y (0.13501144 0.86498856) *
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-50.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-51.png) 

```
##          Prediction
## Reference    N    Y
##         N 3814  127
##         Y  170  693
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   9.381765e-01   7.860827e-01   9.309917e-01   9.448229e-01   8.203580e-01 
## AccuracyPValue  McnemarPValue 
##  2.798570e-127   1.480611e-02
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-52.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-53.png) 

```
##          Prediction
## Reference    N    Y
##         N 1180  318
##         Y   84  146
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   7.673611e-01   2.953321e-01   7.467059e-01   7.871043e-01   8.668981e-01 
## AccuracyPValue  McnemarPValue 
##   1.000000e+00   3.224022e-31 
##                               id                            feats
## 1 Max.cor.Y.rcv.1X1.cp.0###rpart WordCount.root2,NDSSName.my.fctr
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1               0                      0.857                 0.072
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

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-54.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-55.png) 

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
##          NDSSName.my.fctrOpEd#Opnn# NDSSName.my.fctrBsnss#Crsswrds/Gms# 
##                                  55                                  16 
##          NDSSName.my.fctrScnc#Hlth#         NDSSName.my.fctrStyls#U.S.# 
##                                  16                                  12 
## 
## Node number 1: 4804 observations,    complexity param=0.3696408
##   predicted class=N  expected loss=0.179642  P(node) =1
##     class counts:  3941   863
##    probabilities: 0.820 0.180 
##   left son=2 (4367 obs) right son=3 (437 obs)
##   Primary splits:
##       NDSSName.my.fctrOpEd#Opnn#          < 0.5      to the left,  improve=451.59770, (0 missing)
##       NDSSName.my.fctrBsnss#Crsswrds/Gms# < 0.5      to the left,  improve=112.88510, (0 missing)
##       WordCount.root2                     < 25.75849 to the left,  improve=111.17610, (0 missing)
##       NDSSName.my.fctrScnc#Hlth#          < 0.5      to the left,  improve= 99.35206, (0 missing)
##       NDSSName.my.fctrStyls#U.S.#         < 0.5      to the left,  improve= 68.73272, (0 missing)
## 
## Node number 2: 4367 observations,    complexity param=0.09849363
##   predicted class=N  expected loss=0.1110602  P(node) =0.9090341
##     class counts:  3882   485
##    probabilities: 0.889 0.111 
##   left son=4 (4262 obs) right son=5 (105 obs)
##   Primary splits:
##       NDSSName.my.fctrBsnss#Crsswrds/Gms# < 0.5      to the left,  improve=135.55130, (0 missing)
##       NDSSName.my.fctrScnc#Hlth#          < 0.5      to the left,  improve=125.07920, (0 missing)
##       WordCount.root2                     < 25.75849 to the left,  improve= 94.70710, (0 missing)
##       NDSSName.my.fctrStyls#U.S.#         < 0.5      to the left,  improve= 88.56821, (0 missing)
##       NDSSName.my.fctr#Opnn#ThPblcEdtr    < 0.5      to the left,  improve= 18.74400, (0 missing)
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
##       NDSSName.my.fctrScnc#Hlth#       < 0.5      to the left,  improve=132.96710, (0 missing)
##       NDSSName.my.fctrStyls#U.S.#      < 0.5      to the left,  improve= 94.69099, (0 missing)
##       WordCount.root2                  < 26.49528 to the left,  improve= 84.07487, (0 missing)
##       NDSSName.my.fctr#Opnn#ThPblcEdtr < 0.5      to the left,  improve= 19.71762, (0 missing)
##       NDSSName.my.fctrTStyl##          < 0.5      to the right, improve= 10.17000, (0 missing)
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
##       NDSSName.my.fctrStyls#U.S.#      < 0.5      to the left,  improve=102.410700, (0 missing)
##       WordCount.root2                  < 25.01    to the left,  improve= 47.352210, (0 missing)
##       NDSSName.my.fctr#Opnn#ThPblcEdtr < 0.5      to the left,  improve= 20.930810, (0 missing)
##       NDSSName.my.fctrTStyl##          < 0.5      to the right, improve=  5.249425, (0 missing)
##       NDSSName.my.fctrBsnss#Tchnlgy#   < 0.5      to the left,  improve=  2.395935, (0 missing)
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
##    2) NDSSName.my.fctrOpEd#Opnn#< 0.5 4367 485 N (0.88893978 0.11106022)  
##      4) NDSSName.my.fctrBsnss#Crsswrds/Gms#< 0.5 4262 390 N (0.90849366 0.09150634)  
##        8) NDSSName.my.fctrScnc#Hlth#< 0.5 4114 279 N (0.93218279 0.06781721)  
##         16) NDSSName.my.fctrStyls#U.S.#< 0.5 3987 191 N (0.95209431 0.04790569) *
##         17) NDSSName.my.fctrStyls#U.S.#>=0.5 127  39 Y (0.30708661 0.69291339) *
##        9) NDSSName.my.fctrScnc#Hlth#>=0.5 148  37 Y (0.25000000 0.75000000) *
##      5) NDSSName.my.fctrBsnss#Crsswrds/Gms#>=0.5 105  10 Y (0.09523810 0.90476190) *
##    3) NDSSName.my.fctrOpEd#Opnn#>=0.5 437  59 Y (0.13501144 0.86498856) *
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-56.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-57.png) 

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

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-58.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-59.png) 

```
##          Prediction
## Reference    N    Y
##         N 1355  143
##         Y  168   62
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##      0.8200231      0.1825002      0.8010821      0.8378705      0.8668981 
## AccuracyPValue  McnemarPValue 
##      1.0000000      0.1735405 
##                     id                            feats max.nTuningRuns
## 1 Max.cor.Y##rcv#rpart WordCount.root2,NDSSName.my.fctr               5
##   min.elapsedtime.everything min.elapsedtime.final max.AUCpROC.fit
## 1                      2.988                 0.076       0.8709432
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
##                                   label step_major step_minor label_minor
## 5 fit.models_0_Max.cor.Y[rcv.1X1.cp.0|]          1          4       rpart
## 6      fit.models_0_Max.cor.Y.Time.Poly          1          5      glmnet
##       bgn     end elapsed
## 5 169.327 183.585  14.258
## 6 183.586      NA      NA
## [1] "fitting model: Max.cor.Y.Time.Poly##rcv#glmnet"
## [1] "    indep_vars: WordCount.root2,NDSSName.my.fctr,PubDate.day.minutes.poly.1,PubDate.day.minutes.poly.2,PubDate.day.minutes.poly.3,PubDate.day.minutes.poly.4,PubDate.day.minutes.poly.5"
## Aggregating results
## Selecting tuning parameters
## Fitting alpha = 0.55, lambda = 0.0201 on full training set
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-60.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-61.png) 

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
##                         (Intercept)       NDSSName.my.fctr#Opnn#RmFrDbt 
##                        -3.824918286                        -0.730902967 
##    NDSSName.my.fctr#Opnn#ThPblcEdtr NDSSName.my.fctrBsnss#Crsswrds/Gms# 
##                         2.927432404                         3.640472806 
##      NDSSName.my.fctrBsnss#Tchnlgy#    NDSSName.my.fctrFrgn#Wrld#AsPcfc 
##                         0.192915649                        -0.002462006 
##          NDSSName.my.fctrOpEd#Opnn#          NDSSName.my.fctrScnc#Hlth# 
##                         3.940210915                         3.120520277 
##         NDSSName.my.fctrStyls#U.S.#             NDSSName.my.fctrTStyl## 
##                         2.887354523                        -0.349210176 
##          PubDate.day.minutes.poly.1          PubDate.day.minutes.poly.2 
##                        10.171710001                         1.938052563 
##          PubDate.day.minutes.poly.4                     WordCount.root2 
##                         0.422263612                         0.053515500 
## [1] "max lambda < lambdaOpt:"
##                         (Intercept)       NDSSName.my.fctr#Opnn#RmFrDbt 
##                        -3.902999004                        -0.872538727 
##    NDSSName.my.fctr#Opnn#ThPblcEdtr         NDSSName.my.fctr#U.S.#Edctn 
##                         3.038605329                        -0.004334849 
## NDSSName.my.fctrBsnss#Crsswrds/Gms#      NDSSName.my.fctrBsnss#Tchnlgy# 
##                         3.700053147                         0.274702193 
##    NDSSName.my.fctrFrgn#Wrld#AsPcfc          NDSSName.my.fctrOpEd#Opnn# 
##                        -0.061833760                         4.010945376 
##          NDSSName.my.fctrScnc#Hlth#         NDSSName.my.fctrStyls#U.S.# 
##                         3.177975989                         2.950583177 
##             NDSSName.my.fctrTStyl##          PubDate.day.minutes.poly.1 
##                        -0.393389148                        11.082450379 
##          PubDate.day.minutes.poly.2          PubDate.day.minutes.poly.4 
##                         2.849459283                         1.090937007 
##                     WordCount.root2 
##                         0.055710821
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-62.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-63.png) 

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

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-64.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-65.png) 

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
##                                                                                                                                                                     feats
## 1 WordCount.root2,NDSSName.my.fctr,PubDate.day.minutes.poly.1,PubDate.day.minutes.poly.2,PubDate.day.minutes.poly.3,PubDate.day.minutes.poly.4,PubDate.day.minutes.poly.5
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1              25                      4.944                 0.325
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
tune.df = glmnet_tune_models_df,        
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
## 6 fit.models_0_Max.cor.Y.Time.Poly          1          5      glmnet
## 7  fit.models_0_Max.cor.Y.Time.Lag          1          6      glmnet
##       bgn     end elapsed
## 6 183.586 195.343  11.757
## 7 195.344      NA      NA
## [1] "fitting model: Max.cor.Y.Time.Lag##rcv#glmnet"
## [1] "    indep_vars: WordCount.root2,NDSSName.my.fctr,PubDate.last2.log1p,PubDate.last4.log1p,PubDate.last8.log1p,PubDate.last16.log1p,PubDate.last32.log1p"
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

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-66.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-67.png) 

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
##                             (Intercept) 
##                            -3.176434160 
##                 NDSSName.my.fctr#Mltmd# 
##                            -0.126861451 
##           NDSSName.my.fctr#Opnn#RmFrDbt 
##                            -0.486996181 
##        NDSSName.my.fctr#Opnn#ThPblcEdtr 
##                             1.988487294 
##             NDSSName.my.fctr#U.S.#Edctn 
##                            -0.358275828 
##      NDSSName.my.fctrBsnss#BsnssDy#Dlbk 
##                            -0.150520660 
## NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss 
##                            -0.142657391 
##     NDSSName.my.fctrBsnss#Crsswrds/Gms# 
##                             2.404486307 
##              NDSSName.my.fctrCltr#Arts# 
##                            -0.134298496 
##              NDSSName.my.fctrFrgn#Wrld# 
##                            -0.183668961 
##        NDSSName.my.fctrFrgn#Wrld#AsPcfc 
##                            -0.337511604 
##              NDSSName.my.fctrOpEd#Opnn# 
##                             2.519630371 
##              NDSSName.my.fctrScnc#Hlth# 
##                             2.009604889 
##             NDSSName.my.fctrStyls##Fshn 
##                            -0.261634932 
##             NDSSName.my.fctrStyls#U.S.# 
##                             1.848208468 
##                 NDSSName.my.fctrTStyl## 
##                            -0.423871057 
##              NDSSName.my.fctrTrvl#Trvl# 
##                            -0.108409058 
##                     PubDate.last2.log1p 
##                             0.017330994 
##                     PubDate.last4.log1p 
##                             0.030509656 
##                     PubDate.last8.log1p 
##                             0.006803512 
##                         WordCount.root2 
##                             0.031613506 
## [1] "max lambda < lambdaOpt:"
##                             (Intercept) 
##                            -3.298280366 
##                 NDSSName.my.fctr#Mltmd# 
##                            -0.161556150 
##           NDSSName.my.fctr#Opnn#RmFrDbt 
##                            -0.562724348 
##        NDSSName.my.fctr#Opnn#ThPblcEdtr 
##                             2.097804217 
##             NDSSName.my.fctr#U.S.#Edctn 
##                            -0.392971974 
##      NDSSName.my.fctrBsnss#BsnssDy#Dlbk 
##                            -0.160713783 
## NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss 
##                            -0.173587156 
##     NDSSName.my.fctrBsnss#Crsswrds/Gms# 
##                             2.493544311 
##              NDSSName.my.fctrCltr#Arts# 
##                            -0.143142862 
##              NDSSName.my.fctrFrgn#Wrld# 
##                            -0.214101718 
##        NDSSName.my.fctrFrgn#Wrld#AsPcfc 
##                            -0.377780436 
##              NDSSName.my.fctrOpEd#Opnn# 
##                             2.606990045 
##              NDSSName.my.fctrScnc#Hlth# 
##                             2.086724293 
##             NDSSName.my.fctrStyls##Fshn 
##                            -0.300914741 
##             NDSSName.my.fctrStyls#U.S.# 
##                             1.924528272 
##                 NDSSName.my.fctrTStyl## 
##                            -0.450495150 
##              NDSSName.my.fctrTrvl#Trvl# 
##                            -0.140386403 
##                  NDSSName.my.fctrmyOthr 
##                            -0.024631659 
##                     PubDate.last2.log1p 
##                             0.019411174 
##                     PubDate.last4.log1p 
##                             0.033396493 
##                     PubDate.last8.log1p 
##                             0.009744315 
##                         WordCount.root2 
##                             0.033177182
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-68.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-69.png) 

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

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-70.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-71.png) 

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
##                                                                                                                                    feats
## 1 WordCount.root2,NDSSName.my.fctr,PubDate.last2.log1p,PubDate.last4.log1p,PubDate.last8.log1p,PubDate.last16.log1p,PubDate.last32.log1p
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1               5                      3.978                  0.31
##   max.AUCpROC.fit max.Sens.fit max.Spec.fit max.AUCROCR.fit
## 1        0.841928    0.9781781    0.7056779       0.9581739
##   opt.prob.threshold.fit max.f.score.fit max.Accuracy.fit
## 1                    0.3       0.8103347        0.9281861
##   max.AccuracyLower.fit max.AccuracyUpper.fit max.Kappa.fit
## 1              0.925312             0.9396854     0.7362844
##   max.AUCpROC.OOB max.Sens.OOB max.Spec.OOB max.AUCROCR.OOB
## 1       0.5928717    0.9379172    0.2478261       0.8118709
##   opt.prob.threshold.OOB max.f.score.OOB max.Accuracy.OOB
## 1                    0.1       0.4062196         0.646412
##   max.AccuracyLower.OOB max.AccuracyUpper.OOB max.Kappa.OOB
## 1             0.6233476             0.6689781     0.2515036
##   max.AccuracySD.fit max.KappaSD.fit
## 1        0.004492459      0.01977567
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
## 7  fit.models_0_Max.cor.Y.Time.Lag          1          6      glmnet
## 8 fit.models_0_Interact.High.cor.Y          1          7      glmnet
##       bgn     end elapsed
## 7 195.344 206.242  10.899
## 8 206.243      NA      NA
## [1] "fitting model: Interact.High.cor.Y##rcv#glmnet"
## [1] "    indep_vars: WordCount.root2,NDSSName.my.fctr,WordCount.root2:WordCount.root2,WordCount.root2:PubDate.day.minutes.poly.1,WordCount.root2:PubDate.last4.log1p,WordCount.root2:PubDate.last8.log1p,WordCount.root2:PubDate.month.fctr"
## Aggregating results
## Selecting tuning parameters
## Fitting alpha = 0.775, lambda = 0.0201 on full training set
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-72.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-73.png) 

```
##             Length Class      Mode     
## a0            93   -none-     numeric  
## beta        2511   dgCMatrix  S4       
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
## xNames        27   -none-     character
## problemType    1   -none-     character
## tuneValue      2   data.frame list     
## obsLevels      2   -none-     character
## [1] "min lambda > lambdaOpt:"
##                                (Intercept) 
##                               -3.845236612 
##              NDSSName.my.fctr#Opnn#RmFrDbt 
##                               -0.563614974 
##           NDSSName.my.fctr#Opnn#ThPblcEdtr 
##                                2.783870533 
##        NDSSName.my.fctrBsnss#Crsswrds/Gms# 
##                                3.561739787 
##                 NDSSName.my.fctrOpEd#Opnn# 
##                                4.006550818 
##                 NDSSName.my.fctrScnc#Hlth# 
##                                3.128996320 
##                NDSSName.my.fctrStyls#U.S.# 
##                                2.876817882 
##                    NDSSName.my.fctrTStyl## 
##                               -0.111222771 
##                            WordCount.root2 
##                                0.008932084 
## WordCount.root2:PubDate.day.minutes.poly.1 
##                                0.394598722 
##        WordCount.root2:PubDate.last4.log1p 
##                                0.002954057 
##        WordCount.root2:PubDate.last8.log1p 
##                                0.002354682 
## [1] "max lambda < lambdaOpt:"
##                                (Intercept) 
##                               -3.926125697 
##              NDSSName.my.fctr#Opnn#RmFrDbt 
##                               -0.744959921 
##           NDSSName.my.fctr#Opnn#ThPblcEdtr 
##                                2.912506461 
##        NDSSName.my.fctrBsnss#Crsswrds/Gms# 
##                                3.630806510 
##             NDSSName.my.fctrBsnss#Tchnlgy# 
##                                0.084402547 
##                 NDSSName.my.fctrOpEd#Opnn# 
##                                4.076710080 
##                 NDSSName.my.fctrScnc#Hlth# 
##                                3.198327583 
##                NDSSName.my.fctrStyls#U.S.# 
##                                2.943556606 
##                    NDSSName.my.fctrTStyl## 
##                               -0.164715659 
##                            WordCount.root2 
##                                0.009242906 
## WordCount.root2:PubDate.day.minutes.poly.1 
##                                0.450542048 
##        WordCount.root2:PubDate.last4.log1p 
##                                0.003041430 
##        WordCount.root2:PubDate.last8.log1p 
##                                0.002496462
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-74.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-75.png) 

```
##          Prediction
## Reference    N    Y
##         N 3793  148
##         Y  176  687
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   9.325562e-01   7.682400e-01   9.250938e-01   9.394875e-01   8.203580e-01 
## AccuracyPValue  McnemarPValue 
##  9.476102e-114   1.336144e-01
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-76.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-77.png) 

```
##          Prediction
## Reference    N    Y
##         N 1200  298
##         Y   92  138
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   7.743056e-01   2.908255e-01   7.538491e-01   7.938262e-01   8.668981e-01 
## AccuracyPValue  McnemarPValue 
##   1.000000e+00   3.039278e-25 
##                                id
## 1 Interact.High.cor.Y##rcv#glmnet
##                                                                                                                                                                                                                    feats
## 1 WordCount.root2,NDSSName.my.fctr,WordCount.root2:WordCount.root2,WordCount.root2:PubDate.day.minutes.poly.1,WordCount.root2:PubDate.last4.log1p,WordCount.root2:PubDate.last8.log1p,WordCount.root2:PubDate.month.fctr
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1              25                       6.02                 0.408
##   max.AUCpROC.fit max.Sens.fit max.Spec.fit max.AUCROCR.fit
## 1       0.8722119    0.9657447     0.778679         0.95391
##   opt.prob.threshold.fit max.f.score.fit max.Accuracy.fit
## 1                    0.4       0.8091873         0.931585
##   max.AccuracyLower.fit max.AccuracyUpper.fit max.Kappa.fit
## 1             0.9250938             0.9394875     0.7616337
##   max.AUCpROC.OOB max.Sens.OOB max.Spec.OOB max.AUCROCR.OOB
## 1       0.5957392    0.9132176    0.2782609       0.7992251
##   opt.prob.threshold.OOB max.f.score.OOB max.Accuracy.OOB
## 1                    0.1       0.4144144        0.7743056
##   max.AccuracyLower.OOB max.AccuracyUpper.OOB max.Kappa.OOB
## 1             0.7538491             0.7938262     0.2908255
##   max.AccuracySD.fit max.KappaSD.fit
## 1        0.004764279       0.0162102
```

```r
# Low.cor.X
fit.models_0_chunk_df <- myadd_chunk(fit.models_0_chunk_df, 
                        paste0("fit.models_0_", "Low.cor.X"), major.inc = FALSE,
                                     label.minor = "glmnet")
```

```
##                              label step_major step_minor label_minor
## 8 fit.models_0_Interact.High.cor.Y          1          7      glmnet
## 9           fit.models_0_Low.cor.X          1          8      glmnet
##       bgn    end elapsed
## 8 206.243 218.98  12.738
## 9 218.981     NA      NA
```

```r
indep_vars <- subset(glb_feats_df, is.na(cor.high.X) & !nzv & 
                              (exclude.as.feat != 1))[, "id"]  
indep_vars <- myadjust_interaction_feats(indep_vars)
ret_lst <- myfit_mdl(mdl_specs_lst=myinit_mdl_specs_lst(mdl_specs_lst=list(
        id.prefix="Low.cor.X", 
        type=glb_model_type, 
tune.df = glmnet_tune_models_df,        
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
## [1] "    indep_vars: WordCount.root2,NDSSName.my.fctr,PubDate.day.minutes.poly.1,PubDate.wkend,PubDate.day.minutes.poly.4,PubDate.day.minutes.poly.2,PubDate.last4.log1p,PubDate.day.minutes.poly.3,PubDate.last32.log1p,PubDate.month.fctr,.rnorm,PubDate.date.fctr,PubDate.second.fctr,PubDate.minute.fctr,PubDate.wkday.fctr,WordCount.nexp,PubDate.day.minutes.poly.5"
## Aggregating results
## Selecting tuning parameters
## Fitting alpha = 0.1, lambda = 0.0934 on full training set
```

```
## Warning in myfit_mdl(mdl_specs_lst = myinit_mdl_specs_lst(mdl_specs_lst
## = list(id.prefix = "Low.cor.X", : model's bestTune found at an extreme of
## tuneGrid for parameter: alpha
```

```
## Warning in myfit_mdl(mdl_specs_lst = myinit_mdl_specs_lst(mdl_specs_lst
## = list(id.prefix = "Low.cor.X", : model's bestTune found at an extreme of
## tuneGrid for parameter: lambda
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-78.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-79.png) 

```
##             Length Class      Mode     
## a0           100   -none-     numeric  
## beta        5000   dgCMatrix  S4       
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
## xNames        50   -none-     character
## problemType    1   -none-     character
## tuneValue      2   data.frame list     
## obsLevels      2   -none-     character
## [1] "min lambda > lambdaOpt:"
##                             (Intercept) 
##                             -2.99906501 
##                 NDSSName.my.fctr#Mltmd# 
##                             -0.06718324 
##           NDSSName.my.fctr#Opnn#RmFrDbt 
##                             -0.56229342 
##        NDSSName.my.fctr#Opnn#ThPblcEdtr 
##                              1.98953458 
##             NDSSName.my.fctr#U.S.#Edctn 
##                             -0.28336630 
##      NDSSName.my.fctrBsnss#BsnssDy#Dlbk 
##                             -0.14715142 
## NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss 
##                             -0.10505001 
##     NDSSName.my.fctrBsnss#Crsswrds/Gms# 
##                              2.26329365 
##              NDSSName.my.fctrCltr#Arts# 
##                             -0.17558427 
##              NDSSName.my.fctrFrgn#Wrld# 
##                             -0.16747932 
##        NDSSName.my.fctrFrgn#Wrld#AsPcfc 
##                             -0.28174745 
##              NDSSName.my.fctrOpEd#Opnn# 
##                              2.49255592 
##              NDSSName.my.fctrScnc#Hlth# 
##                              2.01908034 
##             NDSSName.my.fctrStyls##Fshn 
##                             -0.23863092 
##             NDSSName.my.fctrStyls#U.S.# 
##                              1.84532884 
##                 NDSSName.my.fctrTStyl## 
##                             -0.43626208 
##              NDSSName.my.fctrTrvl#Trvl# 
##                             -0.12317233 
##              PubDate.day.minutes.poly.1 
##                             10.52089736 
##              PubDate.day.minutes.poly.2 
##                              2.43933073 
##              PubDate.day.minutes.poly.4 
##                              3.62908242 
##                     PubDate.last4.log1p 
##                              0.03007298 
##                           PubDate.wkend 
##                              0.15280167 
##                         WordCount.root2 
##                              0.03121091 
## [1] "max lambda < lambdaOpt:"
##                             (Intercept) 
##                             -3.08485880 
##                 NDSSName.my.fctr#Mltmd# 
##                             -0.09522910 
##           NDSSName.my.fctr#Opnn#RmFrDbt 
##                             -0.64388705 
##        NDSSName.my.fctr#Opnn#ThPblcEdtr 
##                              2.09879679 
##             NDSSName.my.fctr#U.S.#Edctn 
##                             -0.31362838 
##      NDSSName.my.fctrBsnss#BsnssDy#Dlbk 
##                             -0.15761608 
## NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss 
##                             -0.13286636 
##     NDSSName.my.fctrBsnss#Crsswrds/Gms# 
##                              2.34188155 
##              NDSSName.my.fctrCltr#Arts# 
##                             -0.18971959 
##              NDSSName.my.fctrFrgn#Wrld# 
##                             -0.19646484 
##        NDSSName.my.fctrFrgn#Wrld#AsPcfc 
##                             -0.31803524 
##              NDSSName.my.fctrOpEd#Opnn# 
##                              2.57891984 
##              NDSSName.my.fctrScnc#Hlth# 
##                              2.09589193 
##             NDSSName.my.fctrStyls##Fshn 
##                             -0.27680902 
##             NDSSName.my.fctrStyls#U.S.# 
##                              1.92169587 
##                 NDSSName.my.fctrTStyl## 
##                             -0.46355454 
##              NDSSName.my.fctrTrvl#Trvl# 
##                             -0.15595755 
##                  NDSSName.my.fctrmyOthr 
##                             -0.02171951 
##              PubDate.day.minutes.poly.1 
##                             10.96451292 
##              PubDate.day.minutes.poly.2 
##                              2.80569325 
##              PubDate.day.minutes.poly.4 
##                              3.95458868 
##                     PubDate.last4.log1p 
##                              0.03347048 
##                           PubDate.wkend 
##                              0.16348293 
##                         WordCount.root2 
##                              0.03275038
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-80.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-81.png) 

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

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-82.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_0-83.png) 

```
##          Prediction
## Reference   N   Y
##         N 898 600
##         Y  26 204
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   6.377315e-01   2.365595e-01   6.145608e-01   6.604334e-01   8.668981e-01 
## AccuracyPValue  McnemarPValue 
##   1.000000e+00  4.469878e-116 
##                      id
## 1 Low.cor.X##rcv#glmnet
##                                                                                                                                                                                                                                                                                                                                                  feats
## 1 WordCount.root2,NDSSName.my.fctr,PubDate.day.minutes.poly.1,PubDate.wkend,PubDate.day.minutes.poly.4,PubDate.day.minutes.poly.2,PubDate.last4.log1p,PubDate.day.minutes.poly.3,PubDate.last32.log1p,PubDate.month.fctr,.rnorm,PubDate.date.fctr,PubDate.second.fctr,PubDate.minute.fctr,PubDate.wkday.fctr,WordCount.nexp,PubDate.day.minutes.poly.5
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1               5                      5.679                 0.512
##   max.AUCpROC.fit max.Sens.fit max.Spec.fit max.AUCROCR.fit
## 1       0.8383968    0.9769094    0.6998841       0.9582998
##   opt.prob.threshold.fit max.f.score.fit max.Accuracy.fit
## 1                    0.3       0.8084359        0.9261046
##   max.AccuracyLower.fit max.AccuracyUpper.fit max.Kappa.fit
## 1             0.9244394             0.9388938     0.7284544
##   max.AUCpROC.OOB max.Sens.OOB max.Spec.OOB max.AUCROCR.OOB
## 1        0.587685    0.9405874    0.2347826       0.8150026
##   opt.prob.threshold.OOB max.f.score.OOB max.Accuracy.OOB
## 1                    0.1       0.3945841        0.6377315
##   max.AccuracyLower.OOB max.AccuracyUpper.OOB max.Kappa.OOB
## 1             0.6145608             0.6604334     0.2365595
##   max.AccuracySD.fit max.KappaSD.fit
## 1        0.004915343      0.02126499
```

```r
fit.models_0_chunk_df <- 
    myadd_chunk(fit.models_0_chunk_df, "fit.models_0_end", major.inc = FALSE,
                label.minor = "teardown")
```

```
##                     label step_major step_minor label_minor     bgn
## 9  fit.models_0_Low.cor.X          1          8      glmnet 218.981
## 10       fit.models_0_end          1          9    teardown 231.983
##        end elapsed
## 9  231.982  13.001
## 10      NA      NA
```

```r
rm(ret_lst)

glb_chunks_df <- myadd_chunk(glb_chunks_df, "fit.models", major.inc = FALSE)
```

```
##         label step_major step_minor label_minor     bgn     end elapsed
## 10 fit.models          6          0           0  83.208 231.997 148.789
## 11 fit.models          6          1           1 231.997      NA      NA
```


```r
fit.models_1_chunk_df <- myadd_chunk(NULL, "fit.models_1_bgn", label.minor="setup")
```

```
##              label step_major step_minor label_minor     bgn end elapsed
## 1 fit.models_1_bgn          1          0       setup 245.758  NA      NA
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
            tune.df = 
if ((mdl_id_pfx %in% "All.X") && (method %in% "glmnet")) glmnet_tune_models_df else 
                        glb_tune_models_df,
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
## 1   fit.models_1_bgn          1          0       setup 245.758 245.768
## 2 fit.models_1_All.X          1          1       setup 245.769      NA
##   elapsed
## 1   0.011
## 2      NA
##                label step_major step_minor label_minor     bgn     end
## 2 fit.models_1_All.X          1          1       setup 245.769 245.776
## 3 fit.models_1_All.X          1          2      glmnet 245.776      NA
##   elapsed
## 2   0.007
## 3      NA
## [1] "fitting model: All.X##rcv#glmnet"
## [1] "    indep_vars: WordCount.root2,WordCount.log1p,NDSSName.my.fctr,PubDate.day.minutes.poly.1,PubDate.hour.fctr,PubDate.wkend,PubDate.day.minutes.poly.4,PubDate.day.minutes.poly.2,PubDate.last4.log1p,PubDate.last2.log1p,PubDate.last8.log1p,PubDate.last16.log1p,PubDate.day.minutes.poly.3,PubDate.last32.log1p,PubDate.month.fctr,PubDate.juliandate,.rnorm,PubDate.date.fctr,PubDate.second.fctr,PubDate.minute.fctr,PubDate.wkday.fctr,WordCount.nexp,PubDate.day.minutes.poly.5"
## Aggregating results
## Selecting tuning parameters
## Fitting alpha = 0.1, lambda = 0.0934 on full training set
```

```
## Warning in myfit_mdl(mdl_specs_lst = myinit_mdl_specs_lst(mdl_specs_lst
## = list(id.prefix = mdl_id_pfx, : model's bestTune found at an extreme of
## tuneGrid for parameter: alpha
```

```
## Warning in myfit_mdl(mdl_specs_lst = myinit_mdl_specs_lst(mdl_specs_lst
## = list(id.prefix = mdl_id_pfx, : model's bestTune found at an extreme of
## tuneGrid for parameter: lambda
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_1-1.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_1-2.png) 

```
##             Length Class      Mode     
## a0           100   -none-     numeric  
## beta        5700   dgCMatrix  S4       
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
## xNames        57   -none-     character
## problemType    1   -none-     character
## tuneValue      2   data.frame list     
## obsLevels      2   -none-     character
## [1] "min lambda > lambdaOpt:"
##                             (Intercept) 
##                            -3.797789581 
##                 NDSSName.my.fctr#Mltmd# 
##                            -0.033305395 
##           NDSSName.my.fctr#Opnn#RmFrDbt 
##                            -0.609594650 
##        NDSSName.my.fctr#Opnn#ThPblcEdtr 
##                             1.951072023 
##             NDSSName.my.fctr#U.S.#Edctn 
##                            -0.268438616 
##      NDSSName.my.fctrBsnss#BsnssDy#Dlbk 
##                            -0.165907382 
## NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss 
##                            -0.128680248 
##     NDSSName.my.fctrBsnss#Crsswrds/Gms# 
##                             2.226214409 
##              NDSSName.my.fctrCltr#Arts# 
##                            -0.172700226 
##              NDSSName.my.fctrFrgn#Wrld# 
##                            -0.142721283 
##        NDSSName.my.fctrFrgn#Wrld#AsPcfc 
##                            -0.305563165 
##              NDSSName.my.fctrOpEd#Opnn# 
##                             2.487652981 
##              NDSSName.my.fctrScnc#Hlth# 
##                             1.989726118 
##             NDSSName.my.fctrStyls##Fshn 
##                            -0.252459875 
##             NDSSName.my.fctrStyls#U.S.# 
##                             1.824219364 
##                 NDSSName.my.fctrTStyl## 
##                            -0.420594594 
##              NDSSName.my.fctrTrvl#Trvl# 
##                            -0.113478486 
##              PubDate.day.minutes.poly.1 
##                             9.907499536 
##              PubDate.day.minutes.poly.2 
##                             1.591602549 
##              PubDate.day.minutes.poly.4 
##                             4.067701888 
##              PubDate.hour.fctr(15.3,23] 
##                             0.040523622 
##                     PubDate.last2.log1p 
##                             0.012602937 
##                     PubDate.last4.log1p 
##                             0.022795288 
##                     PubDate.last8.log1p 
##                             0.003719588 
##                           PubDate.wkend 
##                             0.143132859 
##                         WordCount.log1p 
##                             0.149525804 
##                         WordCount.root2 
##                             0.023690424 
## [1] "max lambda < lambdaOpt:"
##                             (Intercept) 
##                            -3.930420228 
##                 NDSSName.my.fctr#Mltmd# 
##                            -0.060818787 
##           NDSSName.my.fctr#Opnn#RmFrDbt 
##                            -0.690598861 
##        NDSSName.my.fctr#Opnn#ThPblcEdtr 
##                             2.059371027 
##             NDSSName.my.fctr#U.S.#Edctn 
##                            -0.298341354 
##      NDSSName.my.fctrBsnss#BsnssDy#Dlbk 
##                            -0.177530124 
## NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss 
##                            -0.157260610 
##     NDSSName.my.fctrBsnss#Crsswrds/Gms# 
##                             2.302913242 
##              NDSSName.my.fctrCltr#Arts# 
##                            -0.186887476 
##              NDSSName.my.fctrFrgn#Wrld# 
##                            -0.171078777 
##        NDSSName.my.fctrFrgn#Wrld#AsPcfc 
##                            -0.343168978 
##              NDSSName.my.fctrOpEd#Opnn# 
##                             2.575027796 
##              NDSSName.my.fctrScnc#Hlth# 
##                             2.065753009 
##             NDSSName.my.fctrStyls##Fshn 
##                            -0.291299477 
##             NDSSName.my.fctrStyls#U.S.# 
##                             1.900303104 
##                 NDSSName.my.fctrTStyl## 
##                            -0.447367186 
##              NDSSName.my.fctrTrvl#Trvl# 
##                            -0.146040451 
##                  NDSSName.my.fctrmyOthr 
##                            -0.021887251 
##              PubDate.day.minutes.poly.1 
##                            10.325074080 
##              PubDate.day.minutes.poly.2 
##                             1.878898682 
##              PubDate.day.minutes.poly.4 
##                             4.431162446 
##              PubDate.hour.fctr(15.3,23] 
##                             0.042785355 
##                     PubDate.last2.log1p 
##                             0.014138234 
##                     PubDate.last4.log1p 
##                             0.024615262 
##                     PubDate.last8.log1p 
##                             0.005991056 
##                           PubDate.wkend 
##                             0.151771573 
##                         WordCount.log1p 
##                             0.156106145 
##                         WordCount.root2 
##                             0.024698908
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_1-3.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_1-4.png) 

```
##          Prediction
## Reference    N    Y
##         N 3790  151
##         Y  177  686
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   9.317236e-01   7.655936e-01   9.242213e-01   9.386958e-01   8.203580e-01 
## AccuracyPValue  McnemarPValue 
##  7.764099e-112   1.674653e-01
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_1-5.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.models_1-6.png) 

```
##          Prediction
## Reference   N   Y
##         N 874 624
##         Y  22 208
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   6.261574e-01   2.314268e-01   6.028573e-01   6.490282e-01   8.668981e-01 
## AccuracyPValue  McnemarPValue 
##   1.000000e+00  1.296751e-123 
##                  id
## 1 All.X##rcv#glmnet
##                                                                                                                                                                                                                                                                                                                                                                                                                                                                    feats
## 1 WordCount.root2,WordCount.log1p,NDSSName.my.fctr,PubDate.day.minutes.poly.1,PubDate.hour.fctr,PubDate.wkend,PubDate.day.minutes.poly.4,PubDate.day.minutes.poly.2,PubDate.last4.log1p,PubDate.last2.log1p,PubDate.last8.log1p,PubDate.last16.log1p,PubDate.day.minutes.poly.3,PubDate.last32.log1p,PubDate.month.fctr,PubDate.juliandate,.rnorm,PubDate.date.fctr,PubDate.second.fctr,PubDate.minute.fctr,PubDate.wkday.fctr,WordCount.nexp,PubDate.day.minutes.poly.5
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1               5                      8.046                 0.633
##   max.AUCpROC.fit max.Sens.fit max.Spec.fit max.AUCROCR.fit
## 1       0.8389043    0.9779244    0.6998841       0.9585414
##   opt.prob.threshold.fit max.f.score.fit max.Accuracy.fit
## 1                    0.3       0.8070588        0.9263124
##   max.AccuracyLower.fit max.AccuracyUpper.fit max.Kappa.fit
## 1             0.9242213             0.9386958     0.7277552
##   max.AUCpROC.OOB max.Sens.OOB max.Spec.OOB max.AUCROCR.OOB
## 1        0.587685    0.9405874    0.2347826       0.8157921
##   opt.prob.threshold.OOB max.f.score.OOB max.Accuracy.OOB
## 1                    0.1       0.3917137        0.6261574
##   max.AccuracyLower.OOB max.AccuracyUpper.OOB max.Kappa.OOB
## 1             0.6028573             0.6490282     0.2314268
##   max.AccuracySD.fit max.KappaSD.fit
## 1        0.004624302      0.01929369
```

```r
# Check if other preProcess methods improve model performance
fit.models_1_chunk_df <- 
    myadd_chunk(fit.models_1_chunk_df, "fit.models_1_preProc", major.inc = FALSE,
                label.minor = "preProc")
```

```
##                  label step_major step_minor label_minor     bgn     end
## 3   fit.models_1_All.X          1          2      glmnet 245.776 261.643
## 4 fit.models_1_preProc          1          3     preProc 261.644      NA
##   elapsed
## 3  15.867
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
#                     n_cv_folds=glb_rcv_n_folds, tune_models_df=glb_tune_models_df)
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
##                                                              id
## MFO###myMFO_classfr                         MFO###myMFO_classfr
## Random###myrandom_classfr             Random###myrandom_classfr
## Max.cor.Y.rcv.1X1###glmnet           Max.cor.Y.rcv.1X1###glmnet
## Max.cor.Y.rcv.3X1##rcv#glmnet     Max.cor.Y.rcv.3X1##rcv#glmnet
## Max.cor.Y.rcv.3X3##rcv#glmnet     Max.cor.Y.rcv.3X3##rcv#glmnet
## Max.cor.Y.rcv.3X5##rcv#glmnet     Max.cor.Y.rcv.3X5##rcv#glmnet
## Max.cor.Y.rcv.5X1##rcv#glmnet     Max.cor.Y.rcv.5X1##rcv#glmnet
## Max.cor.Y.rcv.5X3##rcv#glmnet     Max.cor.Y.rcv.5X3##rcv#glmnet
## Max.cor.Y.rcv.5X5##rcv#glmnet     Max.cor.Y.rcv.5X5##rcv#glmnet
## Max.cor.Y.rcv.1X1.cp.0###rpart   Max.cor.Y.rcv.1X1.cp.0###rpart
## Max.cor.Y##rcv#rpart                       Max.cor.Y##rcv#rpart
## Max.cor.Y.Time.Poly##rcv#glmnet Max.cor.Y.Time.Poly##rcv#glmnet
## Max.cor.Y.Time.Lag##rcv#glmnet   Max.cor.Y.Time.Lag##rcv#glmnet
## Interact.High.cor.Y##rcv#glmnet Interact.High.cor.Y##rcv#glmnet
## Low.cor.X##rcv#glmnet                     Low.cor.X##rcv#glmnet
## All.X##rcv#glmnet                             All.X##rcv#glmnet
##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  feats
## MFO###myMFO_classfr                                                                                                                                                                                                                                                                                                                                                                                                                                                                             .rnorm
## Random###myrandom_classfr                                                                                                                                                                                                                                                                                                                                                                                                                                                                       .rnorm
## Max.cor.Y.rcv.1X1###glmnet                                                                                                                                                                                                                                                                                                                                                                                                                                            WordCount.root2,NDSSName.my.fctr
## Max.cor.Y.rcv.3X1##rcv#glmnet                                                                                                                                                                                                                                                                                                                                                                                                                                         WordCount.root2,NDSSName.my.fctr
## Max.cor.Y.rcv.3X3##rcv#glmnet                                                                                                                                                                                                                                                                                                                                                                                                                                         WordCount.root2,NDSSName.my.fctr
## Max.cor.Y.rcv.3X5##rcv#glmnet                                                                                                                                                                                                                                                                                                                                                                                                                                         WordCount.root2,NDSSName.my.fctr
## Max.cor.Y.rcv.5X1##rcv#glmnet                                                                                                                                                                                                                                                                                                                                                                                                                                         WordCount.root2,NDSSName.my.fctr
## Max.cor.Y.rcv.5X3##rcv#glmnet                                                                                                                                                                                                                                                                                                                                                                                                                                         WordCount.root2,NDSSName.my.fctr
## Max.cor.Y.rcv.5X5##rcv#glmnet                                                                                                                                                                                                                                                                                                                                                                                                                                         WordCount.root2,NDSSName.my.fctr
## Max.cor.Y.rcv.1X1.cp.0###rpart                                                                                                                                                                                                                                                                                                                                                                                                                                        WordCount.root2,NDSSName.my.fctr
## Max.cor.Y##rcv#rpart                                                                                                                                                                                                                                                                                                                                                                                                                                                  WordCount.root2,NDSSName.my.fctr
## Max.cor.Y.Time.Poly##rcv#glmnet                                                                                                                                                                                                                                                                                                WordCount.root2,NDSSName.my.fctr,PubDate.day.minutes.poly.1,PubDate.day.minutes.poly.2,PubDate.day.minutes.poly.3,PubDate.day.minutes.poly.4,PubDate.day.minutes.poly.5
## Max.cor.Y.Time.Lag##rcv#glmnet                                                                                                                                                                                                                                                                                                                                  WordCount.root2,NDSSName.my.fctr,PubDate.last2.log1p,PubDate.last4.log1p,PubDate.last8.log1p,PubDate.last16.log1p,PubDate.last32.log1p
## Interact.High.cor.Y##rcv#glmnet                                                                                                                                                                                                                                                 WordCount.root2,NDSSName.my.fctr,WordCount.root2:WordCount.root2,WordCount.root2:PubDate.day.minutes.poly.1,WordCount.root2:PubDate.last4.log1p,WordCount.root2:PubDate.last8.log1p,WordCount.root2:PubDate.month.fctr
## Low.cor.X##rcv#glmnet                                                                                                                             WordCount.root2,NDSSName.my.fctr,PubDate.day.minutes.poly.1,PubDate.wkend,PubDate.day.minutes.poly.4,PubDate.day.minutes.poly.2,PubDate.last4.log1p,PubDate.day.minutes.poly.3,PubDate.last32.log1p,PubDate.month.fctr,.rnorm,PubDate.date.fctr,PubDate.second.fctr,PubDate.minute.fctr,PubDate.wkday.fctr,WordCount.nexp,PubDate.day.minutes.poly.5
## All.X##rcv#glmnet               WordCount.root2,WordCount.log1p,NDSSName.my.fctr,PubDate.day.minutes.poly.1,PubDate.hour.fctr,PubDate.wkend,PubDate.day.minutes.poly.4,PubDate.day.minutes.poly.2,PubDate.last4.log1p,PubDate.last2.log1p,PubDate.last8.log1p,PubDate.last16.log1p,PubDate.day.minutes.poly.3,PubDate.last32.log1p,PubDate.month.fctr,PubDate.juliandate,.rnorm,PubDate.date.fctr,PubDate.second.fctr,PubDate.minute.fctr,PubDate.wkday.fctr,WordCount.nexp,PubDate.day.minutes.poly.5
##                                 max.nTuningRuns min.elapsedtime.everything
## MFO###myMFO_classfr                           0                      0.296
## Random###myrandom_classfr                     0                      0.298
## Max.cor.Y.rcv.1X1###glmnet                    0                      1.024
## Max.cor.Y.rcv.3X1##rcv#glmnet                25                      2.866
## Max.cor.Y.rcv.3X3##rcv#glmnet                25                      4.855
## Max.cor.Y.rcv.3X5##rcv#glmnet                25                      8.311
## Max.cor.Y.rcv.5X1##rcv#glmnet                25                      3.432
## Max.cor.Y.rcv.5X3##rcv#glmnet                25                      6.452
## Max.cor.Y.rcv.5X5##rcv#glmnet                25                      9.150
## Max.cor.Y.rcv.1X1.cp.0###rpart                0                      0.857
## Max.cor.Y##rcv#rpart                          5                      2.988
## Max.cor.Y.Time.Poly##rcv#glmnet              25                      4.944
## Max.cor.Y.Time.Lag##rcv#glmnet                5                      3.978
## Interact.High.cor.Y##rcv#glmnet              25                      6.020
## Low.cor.X##rcv#glmnet                         5                      5.679
## All.X##rcv#glmnet                             5                      8.046
##                                 min.elapsedtime.final max.AUCpROC.fit
## MFO###myMFO_classfr                             0.003       0.5000000
## Random###myrandom_classfr                       0.002       0.4990604
## Max.cor.Y.rcv.1X1###glmnet                      0.276       0.8790544
## Max.cor.Y.rcv.3X1##rcv#glmnet                   0.286       0.8767919
## Max.cor.Y.rcv.3X3##rcv#glmnet                   0.273       0.8767919
## Max.cor.Y.rcv.3X5##rcv#glmnet                   0.277       0.8767919
## Max.cor.Y.rcv.5X1##rcv#glmnet                   0.271       0.8784031
## Max.cor.Y.rcv.5X3##rcv#glmnet                   0.277       0.8784031
## Max.cor.Y.rcv.5X5##rcv#glmnet                   0.269       0.8784031
## Max.cor.Y.rcv.1X1.cp.0###rpart                  0.072       0.8821543
## Max.cor.Y##rcv#rpart                            0.076       0.8709432
## Max.cor.Y.Time.Poly##rcv#glmnet                 0.325       0.8748550
## Max.cor.Y.Time.Lag##rcv#glmnet                  0.310       0.8419280
## Interact.High.cor.Y##rcv#glmnet                 0.408       0.8722119
## Low.cor.X##rcv#glmnet                           0.512       0.8383968
## All.X##rcv#glmnet                               0.633       0.8389043
##                                 max.Sens.fit max.Spec.fit max.AUCROCR.fit
## MFO###myMFO_classfr                1.0000000    0.0000000       0.5000000
## Random###myrandom_classfr          0.8312611    0.1668598       0.4972757
## Max.cor.Y.rcv.1X1###glmnet         0.9632073    0.7949015       0.9608594
## Max.cor.Y.rcv.3X1##rcv#glmnet      0.9644760    0.7891078       0.9582555
## Max.cor.Y.rcv.3X3##rcv#glmnet      0.9644760    0.7891078       0.9582555
## Max.cor.Y.rcv.3X5##rcv#glmnet      0.9644760    0.7891078       0.9582555
## Max.cor.Y.rcv.5X1##rcv#glmnet      0.9642223    0.7925840       0.9607052
## Max.cor.Y.rcv.5X3##rcv#glmnet      0.9642223    0.7925840       0.9607052
## Max.cor.Y.rcv.5X5##rcv#glmnet      0.9642223    0.7925840       0.9607052
## Max.cor.Y.rcv.1X1.cp.0###rpart     0.9705658    0.7937428       0.9504198
## Max.cor.Y##rcv#rpart               0.9632073    0.7786790       0.8746354
## Max.cor.Y.Time.Poly##rcv#glmnet    0.9652372    0.7844728       0.9565422
## Max.cor.Y.Time.Lag##rcv#glmnet     0.9781781    0.7056779       0.9581739
## Interact.High.cor.Y##rcv#glmnet    0.9657447    0.7786790       0.9539100
## Low.cor.X##rcv#glmnet              0.9769094    0.6998841       0.9582998
## All.X##rcv#glmnet                  0.9779244    0.6998841       0.9585414
##                                 opt.prob.threshold.fit max.f.score.fit
## MFO###myMFO_classfr                                0.1       0.3045703
## Random###myrandom_classfr                          0.1       0.3045703
## Max.cor.Y.rcv.1X1###glmnet                         0.5       0.8099174
## Max.cor.Y.rcv.3X1##rcv#glmnet                      0.4       0.8099174
## Max.cor.Y.rcv.3X3##rcv#glmnet                      0.4       0.8099174
## Max.cor.Y.rcv.3X5##rcv#glmnet                      0.4       0.8099174
## Max.cor.Y.rcv.5X1##rcv#glmnet                      0.5       0.8104265
## Max.cor.Y.rcv.5X3##rcv#glmnet                      0.5       0.8104265
## Max.cor.Y.rcv.5X5##rcv#glmnet                      0.5       0.8104265
## Max.cor.Y.rcv.1X1.cp.0###rpart                     0.4       0.8235294
## Max.cor.Y##rcv#rpart                               0.6       0.8000000
## Max.cor.Y.Time.Poly##rcv#glmnet                    0.4       0.8099174
## Max.cor.Y.Time.Lag##rcv#glmnet                     0.3       0.8103347
## Interact.High.cor.Y##rcv#glmnet                    0.4       0.8091873
## Low.cor.X##rcv#glmnet                              0.3       0.8084359
## All.X##rcv#glmnet                                  0.3       0.8070588
##                                 max.Accuracy.fit max.AccuracyLower.fit
## MFO###myMFO_classfr                    0.1796420             0.1688795
## Random###myrandom_classfr              0.1796420             0.1688795
## Max.cor.Y.rcv.1X1###glmnet             0.9329725             0.9255302
## Max.cor.Y.rcv.3X1##rcv#glmnet          0.9335973             0.9255302
## Max.cor.Y.rcv.3X3##rcv#glmnet          0.9333193             0.9255302
## Max.cor.Y.rcv.3X5##rcv#glmnet          0.9332218             0.9255302
## Max.cor.Y.rcv.5X1##rcv#glmnet          0.9331818             0.9259666
## Max.cor.Y.rcv.5X3##rcv#glmnet          0.9333905             0.9259666
## Max.cor.Y.rcv.5X5##rcv#glmnet          0.9331816             0.9259666
## Max.cor.Y.rcv.1X1.cp.0###rpart         0.9381765             0.9309917
## Max.cor.Y##rcv#rpart                   0.9296422             0.9224771
## Max.cor.Y.Time.Poly##rcv#glmnet        0.9322790             0.9255302
## Max.cor.Y.Time.Lag##rcv#glmnet         0.9281861             0.9253120
## Interact.High.cor.Y##rcv#glmnet        0.9315850             0.9250938
## Low.cor.X##rcv#glmnet                  0.9261046             0.9244394
## All.X##rcv#glmnet                      0.9263124             0.9242213
##                                 max.AccuracyUpper.fit max.Kappa.fit
## MFO###myMFO_classfr                         0.1907952     0.0000000
## Random###myrandom_classfr                   0.1907952     0.0000000
## Max.cor.Y.rcv.1X1###glmnet                  0.9398832     0.7692476
## Max.cor.Y.rcv.3X1##rcv#glmnet               0.9398832     0.7691678
## Max.cor.Y.rcv.3X3##rcv#glmnet               0.9398832     0.7690803
## Max.cor.Y.rcv.3X5##rcv#glmnet               0.9398832     0.7686375
## Max.cor.Y.rcv.5X1##rcv#glmnet               0.9402789     0.7689055
## Max.cor.Y.rcv.5X3##rcv#glmnet               0.9402789     0.7698577
## Max.cor.Y.rcv.5X5##rcv#glmnet               0.9402789     0.7691429
## Max.cor.Y.rcv.1X1.cp.0###rpart              0.9448229     0.7860827
## Max.cor.Y##rcv#rpart                        0.9371115     0.7515134
## Max.cor.Y.Time.Poly##rcv#glmnet             0.9398832     0.7643948
## Max.cor.Y.Time.Lag##rcv#glmnet              0.9396854     0.7362844
## Interact.High.cor.Y##rcv#glmnet             0.9394875     0.7616337
## Low.cor.X##rcv#glmnet                       0.9388938     0.7284544
## All.X##rcv#glmnet                           0.9386958     0.7277552
##                                 max.AUCpROC.OOB max.Sens.OOB max.Spec.OOB
## MFO###myMFO_classfr                   0.5000000    1.0000000    0.0000000
## Random###myrandom_classfr             0.5125675    0.8077437    0.2173913
## Max.cor.Y.rcv.1X1###glmnet            0.5962443    0.9098798    0.2826087
## Max.cor.Y.rcv.3X1##rcv#glmnet         0.5962443    0.9098798    0.2826087
## Max.cor.Y.rcv.3X3##rcv#glmnet         0.5962443    0.9098798    0.2826087
## Max.cor.Y.rcv.3X5##rcv#glmnet         0.5962443    0.9098798    0.2826087
## Max.cor.Y.rcv.5X1##rcv#glmnet         0.5962443    0.9098798    0.2826087
## Max.cor.Y.rcv.5X3##rcv#glmnet         0.5962443    0.9098798    0.2826087
## Max.cor.Y.rcv.5X5##rcv#glmnet         0.5962443    0.9098798    0.2826087
## Max.cor.Y.rcv.1X1.cp.0###rpart        0.6174697    0.9218959    0.3130435
## Max.cor.Y##rcv#rpart                  0.5870523    0.9045394    0.2695652
## Max.cor.Y.Time.Poly##rcv#glmnet       0.5965780    0.9105474    0.2826087
## Max.cor.Y.Time.Lag##rcv#glmnet        0.5928717    0.9379172    0.2478261
## Interact.High.cor.Y##rcv#glmnet       0.5957392    0.9132176    0.2782609
## Low.cor.X##rcv#glmnet                 0.5876850    0.9405874    0.2347826
## All.X##rcv#glmnet                     0.5876850    0.9405874    0.2347826
##                                 max.AUCROCR.OOB opt.prob.threshold.OOB
## MFO###myMFO_classfr                   0.5000000                    0.1
## Random###myrandom_classfr             0.4857956                    0.1
## Max.cor.Y.rcv.1X1###glmnet            0.8116126                    0.1
## Max.cor.Y.rcv.3X1##rcv#glmnet         0.8067975                    0.1
## Max.cor.Y.rcv.3X3##rcv#glmnet         0.8067975                    0.1
## Max.cor.Y.rcv.3X5##rcv#glmnet         0.8067975                    0.1
## Max.cor.Y.rcv.5X1##rcv#glmnet         0.8114863                    0.1
## Max.cor.Y.rcv.5X3##rcv#glmnet         0.8114863                    0.1
## Max.cor.Y.rcv.5X5##rcv#glmnet         0.8114863                    0.1
## Max.cor.Y.rcv.1X1.cp.0###rpart        0.7773858                    0.1
## Max.cor.Y##rcv#rpart                  0.5892132                    0.6
## Max.cor.Y.Time.Poly##rcv#glmnet       0.8049472                    0.1
## Max.cor.Y.Time.Lag##rcv#glmnet        0.8118709                    0.1
## Interact.High.cor.Y##rcv#glmnet       0.7992251                    0.1
## Low.cor.X##rcv#glmnet                 0.8150026                    0.1
## All.X##rcv#glmnet                     0.8157921                    0.1
##                                 max.f.score.OOB max.Accuracy.OOB
## MFO###myMFO_classfr                   0.2349336        0.1331019
## Random###myrandom_classfr             0.2349336        0.1331019
## Max.cor.Y.rcv.1X1###glmnet            0.4405405        0.7604167
## Max.cor.Y.rcv.3X1##rcv#glmnet         0.4375839        0.7575231
## Max.cor.Y.rcv.3X3##rcv#glmnet         0.4375839        0.7575231
## Max.cor.Y.rcv.3X5##rcv#glmnet         0.4375839        0.7575231
## Max.cor.Y.rcv.5X1##rcv#glmnet         0.4609375        0.7604167
## Max.cor.Y.rcv.5X3##rcv#glmnet         0.4609375        0.7604167
## Max.cor.Y.rcv.5X5##rcv#glmnet         0.4609375        0.7604167
## Max.cor.Y.rcv.1X1.cp.0###rpart        0.4207493        0.7673611
## Max.cor.Y##rcv#rpart                  0.2850575        0.8200231
## Max.cor.Y.Time.Poly##rcv#glmnet       0.4441261        0.7754630
## Max.cor.Y.Time.Lag##rcv#glmnet        0.4062196        0.6464120
## Interact.High.cor.Y##rcv#glmnet       0.4144144        0.7743056
## Low.cor.X##rcv#glmnet                 0.3945841        0.6377315
## All.X##rcv#glmnet                     0.3917137        0.6261574
##                                 max.AccuracyLower.OOB
## MFO###myMFO_classfr                         0.1174298
## Random###myrandom_classfr                   0.1174298
## Max.cor.Y.rcv.1X1###glmnet                  0.7395703
## Max.cor.Y.rcv.3X1##rcv#glmnet               0.7365992
## Max.cor.Y.rcv.3X3##rcv#glmnet               0.7365992
## Max.cor.Y.rcv.3X5##rcv#glmnet               0.7365992
## Max.cor.Y.rcv.5X1##rcv#glmnet               0.7395703
## Max.cor.Y.rcv.5X3##rcv#glmnet               0.7395703
## Max.cor.Y.rcv.5X5##rcv#glmnet               0.7395703
## Max.cor.Y.rcv.1X1.cp.0###rpart              0.7467059
## Max.cor.Y##rcv#rpart                        0.8010821
## Max.cor.Y.Time.Poly##rcv#glmnet             0.7550404
## Max.cor.Y.Time.Lag##rcv#glmnet              0.6233476
## Interact.High.cor.Y##rcv#glmnet             0.7538491
## Low.cor.X##rcv#glmnet                       0.6145608
## All.X##rcv#glmnet                           0.6028573
##                                 max.AccuracyUpper.OOB max.Kappa.OOB
## MFO###myMFO_classfr                         0.1500310     0.0000000
## Random###myrandom_classfr                   0.1500310     0.0000000
## Max.cor.Y.rcv.1X1###glmnet                  0.7803749     0.3148374
## Max.cor.Y.rcv.3X1##rcv#glmnet               0.7775689     0.3107477
## Max.cor.Y.rcv.3X3##rcv#glmnet               0.7775689     0.3107477
## Max.cor.Y.rcv.3X5##rcv#glmnet               0.7775689     0.3107477
## Max.cor.Y.rcv.5X1##rcv#glmnet               0.7803749     0.3373693
## Max.cor.Y.rcv.5X3##rcv#glmnet               0.7803749     0.3373693
## Max.cor.Y.rcv.5X5##rcv#glmnet               0.7803749     0.3373693
## Max.cor.Y.rcv.1X1.cp.0###rpart              0.7871043     0.2953321
## Max.cor.Y##rcv#rpart                        0.8378705     0.1825002
## Max.cor.Y.Time.Poly##rcv#glmnet             0.7949457     0.3233542
## Max.cor.Y.Time.Lag##rcv#glmnet              0.6689781     0.2515036
## Interact.High.cor.Y##rcv#glmnet             0.7938262     0.2908255
## Low.cor.X##rcv#glmnet                       0.6604334     0.2365595
## All.X##rcv#glmnet                           0.6490282     0.2314268
##                                 max.AccuracySD.fit max.KappaSD.fit
## MFO###myMFO_classfr                             NA              NA
## Random###myrandom_classfr                       NA              NA
## Max.cor.Y.rcv.1X1###glmnet                      NA              NA
## Max.cor.Y.rcv.3X1##rcv#glmnet          0.007015493      0.02403706
## Max.cor.Y.rcv.3X3##rcv#glmnet          0.005178375      0.01754365
## Max.cor.Y.rcv.3X5##rcv#glmnet          0.005396525      0.01835474
## Max.cor.Y.rcv.5X1##rcv#glmnet          0.008837283      0.03133449
## Max.cor.Y.rcv.5X3##rcv#glmnet          0.006138477      0.02161286
## Max.cor.Y.rcv.5X5##rcv#glmnet          0.006213800      0.02210061
## Max.cor.Y.rcv.1X1.cp.0###rpart                  NA              NA
## Max.cor.Y##rcv#rpart                   0.005069520      0.01910910
## Max.cor.Y.Time.Poly##rcv#glmnet        0.005312780      0.01857466
## Max.cor.Y.Time.Lag##rcv#glmnet         0.004492459      0.01977567
## Interact.High.cor.Y##rcv#glmnet        0.004764279      0.01621020
## Low.cor.X##rcv#glmnet                  0.004915343      0.02126499
## All.X##rcv#glmnet                      0.004624302      0.01929369
```

```r
rm(ret_lst)
fit.models_1_chunk_df <- 
    myadd_chunk(fit.models_1_chunk_df, "fit.models_1_end", major.inc = FALSE,
                label.minor = "teardown")
```

```
##                  label step_major step_minor label_minor     bgn     end
## 4 fit.models_1_preProc          1          3     preProc 261.644 261.718
## 5     fit.models_1_end          1          4    teardown 261.719      NA
##   elapsed
## 4   0.074
## 5      NA
```

```r
glb_chunks_df <- myadd_chunk(glb_chunks_df, "fit.models", major.inc = FALSE)
```

```
##         label step_major step_minor label_minor     bgn     end elapsed
## 11 fit.models          6          1           1 231.997 261.728  29.731
## 12 fit.models          6          2           2 261.729      NA      NA
```


```r
fit.models_2_chunk_df <- 
    myadd_chunk(NULL, "fit.models_2_bgn", label.minor = "setup")
```

```
##              label step_major step_minor label_minor     bgn end elapsed
## 1 fit.models_2_bgn          1          0       setup 263.455  NA      NA
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
## Max.cor.Y.rcv.3X1##rcv#glmnet     Max.cor.Y.rcv.3X1##rcv#glmnet
## Max.cor.Y.rcv.3X3##rcv#glmnet     Max.cor.Y.rcv.3X3##rcv#glmnet
## Max.cor.Y.rcv.3X5##rcv#glmnet     Max.cor.Y.rcv.3X5##rcv#glmnet
## Max.cor.Y.rcv.5X1##rcv#glmnet     Max.cor.Y.rcv.5X1##rcv#glmnet
## Max.cor.Y.rcv.5X3##rcv#glmnet     Max.cor.Y.rcv.5X3##rcv#glmnet
## Max.cor.Y.rcv.5X5##rcv#glmnet     Max.cor.Y.rcv.5X5##rcv#glmnet
## Max.cor.Y.rcv.1X1.cp.0###rpart   Max.cor.Y.rcv.1X1.cp.0###rpart
## Max.cor.Y##rcv#rpart                       Max.cor.Y##rcv#rpart
## Max.cor.Y.Time.Poly##rcv#glmnet Max.cor.Y.Time.Poly##rcv#glmnet
## Max.cor.Y.Time.Lag##rcv#glmnet   Max.cor.Y.Time.Lag##rcv#glmnet
## Interact.High.cor.Y##rcv#glmnet Interact.High.cor.Y##rcv#glmnet
## Low.cor.X##rcv#glmnet                     Low.cor.X##rcv#glmnet
## All.X##rcv#glmnet                             All.X##rcv#glmnet
##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  feats
## MFO###myMFO_classfr                                                                                                                                                                                                                                                                                                                                                                                                                                                                             .rnorm
## Random###myrandom_classfr                                                                                                                                                                                                                                                                                                                                                                                                                                                                       .rnorm
## Max.cor.Y.rcv.1X1###glmnet                                                                                                                                                                                                                                                                                                                                                                                                                                            WordCount.root2,NDSSName.my.fctr
## Max.cor.Y.rcv.3X1##rcv#glmnet                                                                                                                                                                                                                                                                                                                                                                                                                                         WordCount.root2,NDSSName.my.fctr
## Max.cor.Y.rcv.3X3##rcv#glmnet                                                                                                                                                                                                                                                                                                                                                                                                                                         WordCount.root2,NDSSName.my.fctr
## Max.cor.Y.rcv.3X5##rcv#glmnet                                                                                                                                                                                                                                                                                                                                                                                                                                         WordCount.root2,NDSSName.my.fctr
## Max.cor.Y.rcv.5X1##rcv#glmnet                                                                                                                                                                                                                                                                                                                                                                                                                                         WordCount.root2,NDSSName.my.fctr
## Max.cor.Y.rcv.5X3##rcv#glmnet                                                                                                                                                                                                                                                                                                                                                                                                                                         WordCount.root2,NDSSName.my.fctr
## Max.cor.Y.rcv.5X5##rcv#glmnet                                                                                                                                                                                                                                                                                                                                                                                                                                         WordCount.root2,NDSSName.my.fctr
## Max.cor.Y.rcv.1X1.cp.0###rpart                                                                                                                                                                                                                                                                                                                                                                                                                                        WordCount.root2,NDSSName.my.fctr
## Max.cor.Y##rcv#rpart                                                                                                                                                                                                                                                                                                                                                                                                                                                  WordCount.root2,NDSSName.my.fctr
## Max.cor.Y.Time.Poly##rcv#glmnet                                                                                                                                                                                                                                                                                                WordCount.root2,NDSSName.my.fctr,PubDate.day.minutes.poly.1,PubDate.day.minutes.poly.2,PubDate.day.minutes.poly.3,PubDate.day.minutes.poly.4,PubDate.day.minutes.poly.5
## Max.cor.Y.Time.Lag##rcv#glmnet                                                                                                                                                                                                                                                                                                                                  WordCount.root2,NDSSName.my.fctr,PubDate.last2.log1p,PubDate.last4.log1p,PubDate.last8.log1p,PubDate.last16.log1p,PubDate.last32.log1p
## Interact.High.cor.Y##rcv#glmnet                                                                                                                                                                                                                                                 WordCount.root2,NDSSName.my.fctr,WordCount.root2:WordCount.root2,WordCount.root2:PubDate.day.minutes.poly.1,WordCount.root2:PubDate.last4.log1p,WordCount.root2:PubDate.last8.log1p,WordCount.root2:PubDate.month.fctr
## Low.cor.X##rcv#glmnet                                                                                                                             WordCount.root2,NDSSName.my.fctr,PubDate.day.minutes.poly.1,PubDate.wkend,PubDate.day.minutes.poly.4,PubDate.day.minutes.poly.2,PubDate.last4.log1p,PubDate.day.minutes.poly.3,PubDate.last32.log1p,PubDate.month.fctr,.rnorm,PubDate.date.fctr,PubDate.second.fctr,PubDate.minute.fctr,PubDate.wkday.fctr,WordCount.nexp,PubDate.day.minutes.poly.5
## All.X##rcv#glmnet               WordCount.root2,WordCount.log1p,NDSSName.my.fctr,PubDate.day.minutes.poly.1,PubDate.hour.fctr,PubDate.wkend,PubDate.day.minutes.poly.4,PubDate.day.minutes.poly.2,PubDate.last4.log1p,PubDate.last2.log1p,PubDate.last8.log1p,PubDate.last16.log1p,PubDate.day.minutes.poly.3,PubDate.last32.log1p,PubDate.month.fctr,PubDate.juliandate,.rnorm,PubDate.date.fctr,PubDate.second.fctr,PubDate.minute.fctr,PubDate.wkday.fctr,WordCount.nexp,PubDate.day.minutes.poly.5
##                                 max.nTuningRuns max.AUCpROC.fit
## MFO###myMFO_classfr                           0       0.5000000
## Random###myrandom_classfr                     0       0.4990604
## Max.cor.Y.rcv.1X1###glmnet                    0       0.8790544
## Max.cor.Y.rcv.3X1##rcv#glmnet                25       0.8767919
## Max.cor.Y.rcv.3X3##rcv#glmnet                25       0.8767919
## Max.cor.Y.rcv.3X5##rcv#glmnet                25       0.8767919
## Max.cor.Y.rcv.5X1##rcv#glmnet                25       0.8784031
## Max.cor.Y.rcv.5X3##rcv#glmnet                25       0.8784031
## Max.cor.Y.rcv.5X5##rcv#glmnet                25       0.8784031
## Max.cor.Y.rcv.1X1.cp.0###rpart                0       0.8821543
## Max.cor.Y##rcv#rpart                          5       0.8709432
## Max.cor.Y.Time.Poly##rcv#glmnet              25       0.8748550
## Max.cor.Y.Time.Lag##rcv#glmnet                5       0.8419280
## Interact.High.cor.Y##rcv#glmnet              25       0.8722119
## Low.cor.X##rcv#glmnet                         5       0.8383968
## All.X##rcv#glmnet                             5       0.8389043
##                                 max.Sens.fit max.Spec.fit max.AUCROCR.fit
## MFO###myMFO_classfr                1.0000000    0.0000000       0.5000000
## Random###myrandom_classfr          0.8312611    0.1668598       0.4972757
## Max.cor.Y.rcv.1X1###glmnet         0.9632073    0.7949015       0.9608594
## Max.cor.Y.rcv.3X1##rcv#glmnet      0.9644760    0.7891078       0.9582555
## Max.cor.Y.rcv.3X3##rcv#glmnet      0.9644760    0.7891078       0.9582555
## Max.cor.Y.rcv.3X5##rcv#glmnet      0.9644760    0.7891078       0.9582555
## Max.cor.Y.rcv.5X1##rcv#glmnet      0.9642223    0.7925840       0.9607052
## Max.cor.Y.rcv.5X3##rcv#glmnet      0.9642223    0.7925840       0.9607052
## Max.cor.Y.rcv.5X5##rcv#glmnet      0.9642223    0.7925840       0.9607052
## Max.cor.Y.rcv.1X1.cp.0###rpart     0.9705658    0.7937428       0.9504198
## Max.cor.Y##rcv#rpart               0.9632073    0.7786790       0.8746354
## Max.cor.Y.Time.Poly##rcv#glmnet    0.9652372    0.7844728       0.9565422
## Max.cor.Y.Time.Lag##rcv#glmnet     0.9781781    0.7056779       0.9581739
## Interact.High.cor.Y##rcv#glmnet    0.9657447    0.7786790       0.9539100
## Low.cor.X##rcv#glmnet              0.9769094    0.6998841       0.9582998
## All.X##rcv#glmnet                  0.9779244    0.6998841       0.9585414
##                                 opt.prob.threshold.fit max.f.score.fit
## MFO###myMFO_classfr                                0.1       0.3045703
## Random###myrandom_classfr                          0.1       0.3045703
## Max.cor.Y.rcv.1X1###glmnet                         0.5       0.8099174
## Max.cor.Y.rcv.3X1##rcv#glmnet                      0.4       0.8099174
## Max.cor.Y.rcv.3X3##rcv#glmnet                      0.4       0.8099174
## Max.cor.Y.rcv.3X5##rcv#glmnet                      0.4       0.8099174
## Max.cor.Y.rcv.5X1##rcv#glmnet                      0.5       0.8104265
## Max.cor.Y.rcv.5X3##rcv#glmnet                      0.5       0.8104265
## Max.cor.Y.rcv.5X5##rcv#glmnet                      0.5       0.8104265
## Max.cor.Y.rcv.1X1.cp.0###rpart                     0.4       0.8235294
## Max.cor.Y##rcv#rpart                               0.6       0.8000000
## Max.cor.Y.Time.Poly##rcv#glmnet                    0.4       0.8099174
## Max.cor.Y.Time.Lag##rcv#glmnet                     0.3       0.8103347
## Interact.High.cor.Y##rcv#glmnet                    0.4       0.8091873
## Low.cor.X##rcv#glmnet                              0.3       0.8084359
## All.X##rcv#glmnet                                  0.3       0.8070588
##                                 max.Accuracy.fit max.Kappa.fit
## MFO###myMFO_classfr                    0.1796420     0.0000000
## Random###myrandom_classfr              0.1796420     0.0000000
## Max.cor.Y.rcv.1X1###glmnet             0.9329725     0.7692476
## Max.cor.Y.rcv.3X1##rcv#glmnet          0.9335973     0.7691678
## Max.cor.Y.rcv.3X3##rcv#glmnet          0.9333193     0.7690803
## Max.cor.Y.rcv.3X5##rcv#glmnet          0.9332218     0.7686375
## Max.cor.Y.rcv.5X1##rcv#glmnet          0.9331818     0.7689055
## Max.cor.Y.rcv.5X3##rcv#glmnet          0.9333905     0.7698577
## Max.cor.Y.rcv.5X5##rcv#glmnet          0.9331816     0.7691429
## Max.cor.Y.rcv.1X1.cp.0###rpart         0.9381765     0.7860827
## Max.cor.Y##rcv#rpart                   0.9296422     0.7515134
## Max.cor.Y.Time.Poly##rcv#glmnet        0.9322790     0.7643948
## Max.cor.Y.Time.Lag##rcv#glmnet         0.9281861     0.7362844
## Interact.High.cor.Y##rcv#glmnet        0.9315850     0.7616337
## Low.cor.X##rcv#glmnet                  0.9261046     0.7284544
## All.X##rcv#glmnet                      0.9263124     0.7277552
##                                 max.AUCpROC.OOB max.Sens.OOB max.Spec.OOB
## MFO###myMFO_classfr                   0.5000000    1.0000000    0.0000000
## Random###myrandom_classfr             0.5125675    0.8077437    0.2173913
## Max.cor.Y.rcv.1X1###glmnet            0.5962443    0.9098798    0.2826087
## Max.cor.Y.rcv.3X1##rcv#glmnet         0.5962443    0.9098798    0.2826087
## Max.cor.Y.rcv.3X3##rcv#glmnet         0.5962443    0.9098798    0.2826087
## Max.cor.Y.rcv.3X5##rcv#glmnet         0.5962443    0.9098798    0.2826087
## Max.cor.Y.rcv.5X1##rcv#glmnet         0.5962443    0.9098798    0.2826087
## Max.cor.Y.rcv.5X3##rcv#glmnet         0.5962443    0.9098798    0.2826087
## Max.cor.Y.rcv.5X5##rcv#glmnet         0.5962443    0.9098798    0.2826087
## Max.cor.Y.rcv.1X1.cp.0###rpart        0.6174697    0.9218959    0.3130435
## Max.cor.Y##rcv#rpart                  0.5870523    0.9045394    0.2695652
## Max.cor.Y.Time.Poly##rcv#glmnet       0.5965780    0.9105474    0.2826087
## Max.cor.Y.Time.Lag##rcv#glmnet        0.5928717    0.9379172    0.2478261
## Interact.High.cor.Y##rcv#glmnet       0.5957392    0.9132176    0.2782609
## Low.cor.X##rcv#glmnet                 0.5876850    0.9405874    0.2347826
## All.X##rcv#glmnet                     0.5876850    0.9405874    0.2347826
##                                 max.AUCROCR.OOB opt.prob.threshold.OOB
## MFO###myMFO_classfr                   0.5000000                    0.1
## Random###myrandom_classfr             0.4857956                    0.1
## Max.cor.Y.rcv.1X1###glmnet            0.8116126                    0.1
## Max.cor.Y.rcv.3X1##rcv#glmnet         0.8067975                    0.1
## Max.cor.Y.rcv.3X3##rcv#glmnet         0.8067975                    0.1
## Max.cor.Y.rcv.3X5##rcv#glmnet         0.8067975                    0.1
## Max.cor.Y.rcv.5X1##rcv#glmnet         0.8114863                    0.1
## Max.cor.Y.rcv.5X3##rcv#glmnet         0.8114863                    0.1
## Max.cor.Y.rcv.5X5##rcv#glmnet         0.8114863                    0.1
## Max.cor.Y.rcv.1X1.cp.0###rpart        0.7773858                    0.1
## Max.cor.Y##rcv#rpart                  0.5892132                    0.6
## Max.cor.Y.Time.Poly##rcv#glmnet       0.8049472                    0.1
## Max.cor.Y.Time.Lag##rcv#glmnet        0.8118709                    0.1
## Interact.High.cor.Y##rcv#glmnet       0.7992251                    0.1
## Low.cor.X##rcv#glmnet                 0.8150026                    0.1
## All.X##rcv#glmnet                     0.8157921                    0.1
##                                 max.f.score.OOB max.Accuracy.OOB
## MFO###myMFO_classfr                   0.2349336        0.1331019
## Random###myrandom_classfr             0.2349336        0.1331019
## Max.cor.Y.rcv.1X1###glmnet            0.4405405        0.7604167
## Max.cor.Y.rcv.3X1##rcv#glmnet         0.4375839        0.7575231
## Max.cor.Y.rcv.3X3##rcv#glmnet         0.4375839        0.7575231
## Max.cor.Y.rcv.3X5##rcv#glmnet         0.4375839        0.7575231
## Max.cor.Y.rcv.5X1##rcv#glmnet         0.4609375        0.7604167
## Max.cor.Y.rcv.5X3##rcv#glmnet         0.4609375        0.7604167
## Max.cor.Y.rcv.5X5##rcv#glmnet         0.4609375        0.7604167
## Max.cor.Y.rcv.1X1.cp.0###rpart        0.4207493        0.7673611
## Max.cor.Y##rcv#rpart                  0.2850575        0.8200231
## Max.cor.Y.Time.Poly##rcv#glmnet       0.4441261        0.7754630
## Max.cor.Y.Time.Lag##rcv#glmnet        0.4062196        0.6464120
## Interact.High.cor.Y##rcv#glmnet       0.4144144        0.7743056
## Low.cor.X##rcv#glmnet                 0.3945841        0.6377315
## All.X##rcv#glmnet                     0.3917137        0.6261574
##                                 max.Kappa.OOB inv.elapsedtime.everything
## MFO###myMFO_classfr                 0.0000000                  3.3783784
## Random###myrandom_classfr           0.0000000                  3.3557047
## Max.cor.Y.rcv.1X1###glmnet          0.3148374                  0.9765625
## Max.cor.Y.rcv.3X1##rcv#glmnet       0.3107477                  0.3489184
## Max.cor.Y.rcv.3X3##rcv#glmnet       0.3107477                  0.2059732
## Max.cor.Y.rcv.3X5##rcv#glmnet       0.3107477                  0.1203225
## Max.cor.Y.rcv.5X1##rcv#glmnet       0.3373693                  0.2913753
## Max.cor.Y.rcv.5X3##rcv#glmnet       0.3373693                  0.1549907
## Max.cor.Y.rcv.5X5##rcv#glmnet       0.3373693                  0.1092896
## Max.cor.Y.rcv.1X1.cp.0###rpart      0.2953321                  1.1668611
## Max.cor.Y##rcv#rpart                0.1825002                  0.3346720
## Max.cor.Y.Time.Poly##rcv#glmnet     0.3233542                  0.2022654
## Max.cor.Y.Time.Lag##rcv#glmnet      0.2515036                  0.2513826
## Interact.High.cor.Y##rcv#glmnet     0.2908255                  0.1661130
## Low.cor.X##rcv#glmnet               0.2365595                  0.1760873
## All.X##rcv#glmnet                   0.2314268                  0.1242854
##                                 inv.elapsedtime.final
## MFO###myMFO_classfr                        333.333333
## Random###myrandom_classfr                  500.000000
## Max.cor.Y.rcv.1X1###glmnet                   3.623188
## Max.cor.Y.rcv.3X1##rcv#glmnet                3.496503
## Max.cor.Y.rcv.3X3##rcv#glmnet                3.663004
## Max.cor.Y.rcv.3X5##rcv#glmnet                3.610108
## Max.cor.Y.rcv.5X1##rcv#glmnet                3.690037
## Max.cor.Y.rcv.5X3##rcv#glmnet                3.610108
## Max.cor.Y.rcv.5X5##rcv#glmnet                3.717472
## Max.cor.Y.rcv.1X1.cp.0###rpart              13.888889
## Max.cor.Y##rcv#rpart                        13.157895
## Max.cor.Y.Time.Poly##rcv#glmnet              3.076923
## Max.cor.Y.Time.Lag##rcv#glmnet               3.225806
## Interact.High.cor.Y##rcv#glmnet              2.450980
## Low.cor.X##rcv#glmnet                        1.953125
## All.X##rcv#glmnet                            1.579779
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
## 16. Consider specifying shapes manually if you must have them.
```

```
## Warning: Removed 200 rows containing missing values (geom_point).
```

```
## Warning in RColorBrewer::brewer.pal(n, pal): n too large, allowed maximum for palette Set1 is 9
## Returning the palette you asked for with that many colors
```

```
## Warning: The shape palette can deal with a maximum of 6 discrete values
## because more than 6 becomes difficult to discriminate; you have
## 16. Consider specifying shapes manually if you must have them.
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
## Interact.High.cor.Y##rcv#glmnet Interact.High.cor.Y##rcv#glmnet
## Max.cor.Y.rcv.1X1.cp.0###rpart   Max.cor.Y.rcv.1X1.cp.0###rpart
## Max.cor.Y.rcv.1X1###glmnet           Max.cor.Y.rcv.1X1###glmnet
## Max.cor.Y.rcv.5X3##rcv#glmnet     Max.cor.Y.rcv.5X3##rcv#glmnet
## Max.cor.Y.rcv.5X1##rcv#glmnet     Max.cor.Y.rcv.5X1##rcv#glmnet
## Max.cor.Y.rcv.5X5##rcv#glmnet     Max.cor.Y.rcv.5X5##rcv#glmnet
## Max.cor.Y.rcv.3X1##rcv#glmnet     Max.cor.Y.rcv.3X1##rcv#glmnet
## Max.cor.Y.rcv.3X3##rcv#glmnet     Max.cor.Y.rcv.3X3##rcv#glmnet
## Max.cor.Y.rcv.3X5##rcv#glmnet     Max.cor.Y.rcv.3X5##rcv#glmnet
## Max.cor.Y.Time.Lag##rcv#glmnet   Max.cor.Y.Time.Lag##rcv#glmnet
## Low.cor.X##rcv#glmnet                     Low.cor.X##rcv#glmnet
## All.X##rcv#glmnet                             All.X##rcv#glmnet
## MFO###myMFO_classfr                         MFO###myMFO_classfr
## Random###myrandom_classfr             Random###myrandom_classfr
##                                 max.Accuracy.OOB max.AUCROCR.OOB
## Max.cor.Y##rcv#rpart                   0.8200231       0.5892132
## Max.cor.Y.Time.Poly##rcv#glmnet        0.7754630       0.8049472
## Interact.High.cor.Y##rcv#glmnet        0.7743056       0.7992251
## Max.cor.Y.rcv.1X1.cp.0###rpart         0.7673611       0.7773858
## Max.cor.Y.rcv.1X1###glmnet             0.7604167       0.8116126
## Max.cor.Y.rcv.5X3##rcv#glmnet          0.7604167       0.8114863
## Max.cor.Y.rcv.5X1##rcv#glmnet          0.7604167       0.8114863
## Max.cor.Y.rcv.5X5##rcv#glmnet          0.7604167       0.8114863
## Max.cor.Y.rcv.3X1##rcv#glmnet          0.7575231       0.8067975
## Max.cor.Y.rcv.3X3##rcv#glmnet          0.7575231       0.8067975
## Max.cor.Y.rcv.3X5##rcv#glmnet          0.7575231       0.8067975
## Max.cor.Y.Time.Lag##rcv#glmnet         0.6464120       0.8118709
## Low.cor.X##rcv#glmnet                  0.6377315       0.8150026
## All.X##rcv#glmnet                      0.6261574       0.8157921
## MFO###myMFO_classfr                    0.1331019       0.5000000
## Random###myrandom_classfr              0.1331019       0.4857956
##                                 max.AUCpROC.OOB max.Accuracy.fit
## Max.cor.Y##rcv#rpart                  0.5870523        0.9296422
## Max.cor.Y.Time.Poly##rcv#glmnet       0.5965780        0.9322790
## Interact.High.cor.Y##rcv#glmnet       0.5957392        0.9315850
## Max.cor.Y.rcv.1X1.cp.0###rpart        0.6174697        0.9381765
## Max.cor.Y.rcv.1X1###glmnet            0.5962443        0.9329725
## Max.cor.Y.rcv.5X3##rcv#glmnet         0.5962443        0.9333905
## Max.cor.Y.rcv.5X1##rcv#glmnet         0.5962443        0.9331818
## Max.cor.Y.rcv.5X5##rcv#glmnet         0.5962443        0.9331816
## Max.cor.Y.rcv.3X1##rcv#glmnet         0.5962443        0.9335973
## Max.cor.Y.rcv.3X3##rcv#glmnet         0.5962443        0.9333193
## Max.cor.Y.rcv.3X5##rcv#glmnet         0.5962443        0.9332218
## Max.cor.Y.Time.Lag##rcv#glmnet        0.5928717        0.9281861
## Low.cor.X##rcv#glmnet                 0.5876850        0.9261046
## All.X##rcv#glmnet                     0.5876850        0.9263124
## MFO###myMFO_classfr                   0.5000000        0.1796420
## Random###myrandom_classfr             0.5125675        0.1796420
##                                 opt.prob.threshold.fit
## Max.cor.Y##rcv#rpart                               0.6
## Max.cor.Y.Time.Poly##rcv#glmnet                    0.4
## Interact.High.cor.Y##rcv#glmnet                    0.4
## Max.cor.Y.rcv.1X1.cp.0###rpart                     0.4
## Max.cor.Y.rcv.1X1###glmnet                         0.5
## Max.cor.Y.rcv.5X3##rcv#glmnet                      0.5
## Max.cor.Y.rcv.5X1##rcv#glmnet                      0.5
## Max.cor.Y.rcv.5X5##rcv#glmnet                      0.5
## Max.cor.Y.rcv.3X1##rcv#glmnet                      0.4
## Max.cor.Y.rcv.3X3##rcv#glmnet                      0.4
## Max.cor.Y.rcv.3X5##rcv#glmnet                      0.4
## Max.cor.Y.Time.Lag##rcv#glmnet                     0.3
## Low.cor.X##rcv#glmnet                              0.3
## All.X##rcv#glmnet                                  0.3
## MFO###myMFO_classfr                                0.1
## Random###myrandom_classfr                          0.1
##                                 opt.prob.threshold.OOB
## Max.cor.Y##rcv#rpart                               0.6
## Max.cor.Y.Time.Poly##rcv#glmnet                    0.1
## Interact.High.cor.Y##rcv#glmnet                    0.1
## Max.cor.Y.rcv.1X1.cp.0###rpart                     0.1
## Max.cor.Y.rcv.1X1###glmnet                         0.1
## Max.cor.Y.rcv.5X3##rcv#glmnet                      0.1
## Max.cor.Y.rcv.5X1##rcv#glmnet                      0.1
## Max.cor.Y.rcv.5X5##rcv#glmnet                      0.1
## Max.cor.Y.rcv.3X1##rcv#glmnet                      0.1
## Max.cor.Y.rcv.3X3##rcv#glmnet                      0.1
## Max.cor.Y.rcv.3X5##rcv#glmnet                      0.1
## Max.cor.Y.Time.Lag##rcv#glmnet                     0.1
## Low.cor.X##rcv#glmnet                              0.1
## All.X##rcv#glmnet                                  0.1
## MFO###myMFO_classfr                                0.1
## Random###myrandom_classfr                          0.1
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
## 16. Consider specifying shapes manually if you must have them.
```

```
## Warning: Removed 70 rows containing missing values (geom_point).
```

```
## Warning in RColorBrewer::brewer.pal(n, pal): n too large, allowed maximum for palette Set1 is 9
## Returning the palette you asked for with that many colors
```

```
## Warning: The shape palette can deal with a maximum of 6 discrete values
## because more than 6 becomes difficult to discriminate; you have
## 16. Consider specifying shapes manually if you must have them.
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
## <environment: 0x7f8ff66ebce0>
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
## a0           100   -none-     numeric  
## beta        5700   dgCMatrix  S4       
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
## xNames        57   -none-     character
## problemType    1   -none-     character
## tuneValue      2   data.frame list     
## obsLevels      2   -none-     character
## [1] "min lambda > lambdaOpt:"
##                             (Intercept) 
##                            -3.797789581 
##                 NDSSName.my.fctr#Mltmd# 
##                            -0.033305395 
##           NDSSName.my.fctr#Opnn#RmFrDbt 
##                            -0.609594650 
##        NDSSName.my.fctr#Opnn#ThPblcEdtr 
##                             1.951072023 
##             NDSSName.my.fctr#U.S.#Edctn 
##                            -0.268438616 
##      NDSSName.my.fctrBsnss#BsnssDy#Dlbk 
##                            -0.165907382 
## NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss 
##                            -0.128680248 
##     NDSSName.my.fctrBsnss#Crsswrds/Gms# 
##                             2.226214409 
##              NDSSName.my.fctrCltr#Arts# 
##                            -0.172700226 
##              NDSSName.my.fctrFrgn#Wrld# 
##                            -0.142721283 
##        NDSSName.my.fctrFrgn#Wrld#AsPcfc 
##                            -0.305563165 
##              NDSSName.my.fctrOpEd#Opnn# 
##                             2.487652981 
##              NDSSName.my.fctrScnc#Hlth# 
##                             1.989726118 
##             NDSSName.my.fctrStyls##Fshn 
##                            -0.252459875 
##             NDSSName.my.fctrStyls#U.S.# 
##                             1.824219364 
##                 NDSSName.my.fctrTStyl## 
##                            -0.420594594 
##              NDSSName.my.fctrTrvl#Trvl# 
##                            -0.113478486 
##              PubDate.day.minutes.poly.1 
##                             9.907499536 
##              PubDate.day.minutes.poly.2 
##                             1.591602549 
##              PubDate.day.minutes.poly.4 
##                             4.067701888 
##              PubDate.hour.fctr(15.3,23] 
##                             0.040523622 
##                     PubDate.last2.log1p 
##                             0.012602937 
##                     PubDate.last4.log1p 
##                             0.022795288 
##                     PubDate.last8.log1p 
##                             0.003719588 
##                           PubDate.wkend 
##                             0.143132859 
##                         WordCount.log1p 
##                             0.149525804 
##                         WordCount.root2 
##                             0.023690424 
## [1] "max lambda < lambdaOpt:"
##                             (Intercept) 
##                            -3.930420228 
##                 NDSSName.my.fctr#Mltmd# 
##                            -0.060818787 
##           NDSSName.my.fctr#Opnn#RmFrDbt 
##                            -0.690598861 
##        NDSSName.my.fctr#Opnn#ThPblcEdtr 
##                             2.059371027 
##             NDSSName.my.fctr#U.S.#Edctn 
##                            -0.298341354 
##      NDSSName.my.fctrBsnss#BsnssDy#Dlbk 
##                            -0.177530124 
## NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss 
##                            -0.157260610 
##     NDSSName.my.fctrBsnss#Crsswrds/Gms# 
##                             2.302913242 
##              NDSSName.my.fctrCltr#Arts# 
##                            -0.186887476 
##              NDSSName.my.fctrFrgn#Wrld# 
##                            -0.171078777 
##        NDSSName.my.fctrFrgn#Wrld#AsPcfc 
##                            -0.343168978 
##              NDSSName.my.fctrOpEd#Opnn# 
##                             2.575027796 
##              NDSSName.my.fctrScnc#Hlth# 
##                             2.065753009 
##             NDSSName.my.fctrStyls##Fshn 
##                            -0.291299477 
##             NDSSName.my.fctrStyls#U.S.# 
##                             1.900303104 
##                 NDSSName.my.fctrTStyl## 
##                            -0.447367186 
##              NDSSName.my.fctrTrvl#Trvl# 
##                            -0.146040451 
##                  NDSSName.my.fctrmyOthr 
##                            -0.021887251 
##              PubDate.day.minutes.poly.1 
##                            10.325074080 
##              PubDate.day.minutes.poly.2 
##                             1.878898682 
##              PubDate.day.minutes.poly.4 
##                             4.431162446 
##              PubDate.hour.fctr(15.3,23] 
##                             0.042785355 
##                     PubDate.last2.log1p 
##                             0.014138234 
##                     PubDate.last4.log1p 
##                             0.024615262 
##                     PubDate.last8.log1p 
##                             0.005991056 
##                           PubDate.wkend 
##                             0.151771573 
##                         WordCount.log1p 
##                             0.156106145 
##                         WordCount.root2 
##                             0.024698908
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
##                                                imp All.X##rcv#glmnet.imp
## PubDate.day.minutes.poly.1              100.000000            100.000000
## PubDate.day.minutes.poly.4               46.118345             46.118345
## NDSSName.my.fctrOpEd#Opnn#               29.608813             29.608813
## NDSSName.my.fctrBsnss#Crsswrds/Gms#      27.135657             27.135657
## NDSSName.my.fctrScnc#Hlth#               24.964890             24.964890
## NDSSName.my.fctr#Opnn#ThPblcEdtr         24.849247             24.849247
## NDSSName.my.fctrStyls#U.S.#              23.449567             23.449567
## PubDate.day.minutes.poly.2               22.879214             22.879214
## WordCount.log1p                           7.599061              7.599061
## PubDate.wkend                             7.555716              7.555716
## PubDate.hour.fctr(15.3,23]                6.568901              6.568901
## WordCount.root2                           6.405483              6.405483
## PubDate.last4.log1p                       6.403279              6.403279
## PubDate.last2.log1p                       6.307833              6.307833
## PubDate.last8.log1p                       6.231915              6.231915
## .rnorm                                    6.181073              6.181073
## NDSSName.my.fctrBsnss#Tchnlgy#            6.181073              6.181073
## NDSSName.my.fctrCltr##                    6.181073              6.181073
## NDSSName.my.fctrMtr#N.Y./Rgn#             6.181073              6.181073
## PubDate.date.fctr(7,13]                   6.181073              6.181073
## PubDate.date.fctr(13,19]                  6.181073              6.181073
## PubDate.date.fctr(19,25]                  6.181073              6.181073
## PubDate.date.fctr(25,31]                  6.181073              6.181073
## PubDate.day.minutes.poly.3                6.181073              6.181073
## PubDate.day.minutes.poly.5                6.181073              6.181073
## PubDate.hour.fctr(7.67,15.3]              6.181073              6.181073
## PubDate.juliandate                        6.181073              6.181073
## PubDate.last16.log1p                      6.181073              6.181073
## PubDate.last32.log1p                      6.181073              6.181073
## PubDate.minute.fctr(14.8,29.5]            6.181073              6.181073
## PubDate.minute.fctr(29.5,44.2]            6.181073              6.181073
## PubDate.minute.fctr(44.2,59.1]            6.181073              6.181073
## PubDate.month.fctr10                      6.181073              6.181073
## PubDate.month.fctr11                      6.181073              6.181073
## PubDate.month.fctr12                      6.181073              6.181073
## PubDate.second.fctr(14.8,29.5]            6.181073              6.181073
## PubDate.second.fctr(29.5,44.2]            6.181073              6.181073
## PubDate.second.fctr(44.2,59.1]            6.181073              6.181073
## PubDate.wkday.fctr1                       6.181073              6.181073
## PubDate.wkday.fctr2                       6.181073              6.181073
## PubDate.wkday.fctr3                       6.181073              6.181073
## PubDate.wkday.fctr4                       6.181073              6.181073
## PubDate.wkday.fctr5                       6.181073              6.181073
## PubDate.wkday.fctr6                       6.181073              6.181073
## WordCount.nexp                            6.181073              6.181073
## NDSSName.my.fctrmyOthr                    6.019416              6.019416
## NDSSName.my.fctr#Mltmd#                   5.672845              5.672845
## NDSSName.my.fctrTrvl#Trvl#                4.901316              4.901316
## NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss   4.791503              4.791503
## NDSSName.my.fctrFrgn#Wrld#                4.664559              4.664559
## NDSSName.my.fctrBsnss#BsnssDy#Dlbk        4.575817              4.575817
## NDSSName.my.fctrCltr#Arts#                4.494666              4.494666
## NDSSName.my.fctrStyls##Fshn               3.582132              3.582132
## NDSSName.my.fctr#U.S.#Edctn               3.501802              3.501802
## NDSSName.my.fctrFrgn#Wrld#AsPcfc          3.104915              3.104915
## NDSSName.my.fctrTStyl##                   2.131448              2.131448
## NDSSName.my.fctr#Opnn#RmFrDbt             0.000000              0.000000
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
## 1     3918         N                       0.03938738
## 2     2555         N                       0.02746354
## 3      302         N                       0.24364717
##   Pplr.fctr.All.X..rcv.glmnet Pplr.fctr.All.X..rcv.glmnet.err
## 1                           N                           FALSE
## 2                           N                           FALSE
## 3                           Y                            TRUE
##   Pplr.fctr.All.X..rcv.glmnet.err.abs Pplr.fctr.All.X..rcv.glmnet.is.acc
## 1                          0.03938738                               TRUE
## 2                          0.02746354                               TRUE
## 3                          0.24364717                              FALSE
##   Pplr.fctr.All.X..rcv.glmnet.accurate Pplr.fctr.All.X..rcv.glmnet.error
## 1                                 TRUE                         0.0000000
## 2                                 TRUE                         0.0000000
## 3                                FALSE                         0.1436472
##   .label
## 1   3918
## 2   2555
## 3    302
## [1] "Inaccurate: "
##   UniqueID Pplr.fctr Pplr.fctr.All.X..rcv.glmnet.prob
## 1      172         Y                       0.06329231
## 2     3554         Y                       0.06519608
## 3       92         Y                       0.06728628
## 4     3076         Y                       0.07036009
## 5     6354         Y                       0.07093837
## 6     4020         Y                       0.07165963
##   Pplr.fctr.All.X..rcv.glmnet Pplr.fctr.All.X..rcv.glmnet.err
## 1                           N                            TRUE
## 2                           N                            TRUE
## 3                           N                            TRUE
## 4                           N                            TRUE
## 5                           N                            TRUE
## 6                           N                            TRUE
##   Pplr.fctr.All.X..rcv.glmnet.err.abs Pplr.fctr.All.X..rcv.glmnet.is.acc
## 1                           0.9367077                              FALSE
## 2                           0.9348039                              FALSE
## 3                           0.9327137                              FALSE
## 4                           0.9296399                              FALSE
## 5                           0.9290616                              FALSE
## 6                           0.9283404                              FALSE
##   Pplr.fctr.All.X..rcv.glmnet.accurate Pplr.fctr.All.X..rcv.glmnet.error
## 1                                FALSE                       -0.03670769
## 2                                FALSE                       -0.03480392
## 3                                FALSE                       -0.03271372
## 4                                FALSE                       -0.02963991
## 5                                FALSE                       -0.02906163
## 6                                FALSE                       -0.02834037
##     UniqueID Pplr.fctr Pplr.fctr.All.X..rcv.glmnet.prob
## 156     1774         N                        0.1120934
## 332     4664         N                        0.1360889
## 361     6480         N                        0.1410233
## 396      508         N                        0.1509912
## 400     1922         N                        0.1520510
## 635      483         N                        0.7423412
##     Pplr.fctr.All.X..rcv.glmnet Pplr.fctr.All.X..rcv.glmnet.err
## 156                           Y                            TRUE
## 332                           Y                            TRUE
## 361                           Y                            TRUE
## 396                           Y                            TRUE
## 400                           Y                            TRUE
## 635                           Y                            TRUE
##     Pplr.fctr.All.X..rcv.glmnet.err.abs Pplr.fctr.All.X..rcv.glmnet.is.acc
## 156                           0.1120934                              FALSE
## 332                           0.1360889                              FALSE
## 361                           0.1410233                              FALSE
## 396                           0.1509912                              FALSE
## 400                           0.1520510                              FALSE
## 635                           0.7423412                              FALSE
##     Pplr.fctr.All.X..rcv.glmnet.accurate Pplr.fctr.All.X..rcv.glmnet.error
## 156                                FALSE                        0.01209336
## 332                                FALSE                        0.03608890
## 361                                FALSE                        0.04102331
## 396                                FALSE                        0.05099123
## 400                                FALSE                        0.05205101
## 635                                FALSE                        0.64234122
##     UniqueID Pplr.fctr Pplr.fctr.All.X..rcv.glmnet.prob
## 641      770         N                        0.7842764
## 642      221         N                        0.7905685
## 643      472         N                        0.7907410
## 644     1448         N                        0.8011158
## 645     3590         N                        0.8019486
## 646     2995         N                        0.8069664
##     Pplr.fctr.All.X..rcv.glmnet Pplr.fctr.All.X..rcv.glmnet.err
## 641                           Y                            TRUE
## 642                           Y                            TRUE
## 643                           Y                            TRUE
## 644                           Y                            TRUE
## 645                           Y                            TRUE
## 646                           Y                            TRUE
##     Pplr.fctr.All.X..rcv.glmnet.err.abs Pplr.fctr.All.X..rcv.glmnet.is.acc
## 641                           0.7842764                              FALSE
## 642                           0.7905685                              FALSE
## 643                           0.7907410                              FALSE
## 644                           0.8011158                              FALSE
## 645                           0.8019486                              FALSE
## 646                           0.8069664                              FALSE
##     Pplr.fctr.All.X..rcv.glmnet.accurate Pplr.fctr.All.X..rcv.glmnet.error
## 641                                FALSE                         0.6842764
## 642                                FALSE                         0.6905685
## 643                                FALSE                         0.6907410
## 644                                FALSE                         0.7011158
## 645                                FALSE                         0.7019486
## 646                                FALSE                         0.7069664
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
##                                NDSSName.my.fctr .n.OOB .n.Fit .n.Tst
## OpEd#Opnn#                           OpEd#Opnn#     89    437    164
## #Opnn#ThPblcEdtr               #Opnn#ThPblcEdtr      4     16     10
## Styls#U.S.#                         Styls#U.S.#     50    127     61
## Bsnss#Crsswrds/Gms#         Bsnss#Crsswrds/Gms#     18    105     42
## Scnc#Hlth#                           Scnc#Hlth#     48    148     57
## Bsnss#Tchnlgy#                   Bsnss#Tchnlgy#    126    213    114
## ##                                           ##    371    913    342
## Bsnss#BsnssDy#Dlbk           Bsnss#BsnssDy#Dlbk    323    629    304
## Mtr#N.Y./Rgn#                     Mtr#N.Y./Rgn#     70    128     67
## Cltr#Arts#                           Cltr#Arts#    185    490    174
## #Opnn#RmFrDbt                     #Opnn#RmFrDbt     20     42     20
## Styls##Fshn                         Styls##Fshn     15    104     15
## Bsnss#BsnssDy#SmllBsnss Bsnss#BsnssDy#SmllBsnss     40    100     41
## myOthr                                   myOthr      5     33      5
## Trvl#Trvl#                           Trvl#Trvl#     34     83     35
## Cltr##                                   Cltr##      1     NA     70
## Frgn#Wrld#AsPcfc               Frgn#Wrld#AsPcfc     53    150     56
## #Mltmd#                                 #Mltmd#     49     92     52
## TStyl##                                 TStyl##    101    623    105
## #U.S.#Edctn                         #U.S.#Edctn     82    243     89
## Frgn#Wrld#                           Frgn#Wrld#     44    128     47
##                         .freqRatio.Fit .freqRatio.OOB .freqRatio.Tst
## OpEd#Opnn#                 0.090965862   0.0515046296    0.087700535
## #Opnn#ThPblcEdtr           0.003330558   0.0023148148    0.005347594
## Styls#U.S.#                0.026436303   0.0289351852    0.032620321
## Bsnss#Crsswrds/Gms#        0.021856786   0.0104166667    0.022459893
## Scnc#Hlth#                 0.030807660   0.0277777778    0.030481283
## Bsnss#Tchnlgy#             0.044338052   0.0729166667    0.060962567
## ##                         0.190049958   0.2146990741    0.182887701
## Bsnss#BsnssDy#Dlbk         0.130932556   0.1869212963    0.162566845
## Mtr#N.Y./Rgn#              0.026644463   0.0405092593    0.035828877
## Cltr#Arts#                 0.101998335   0.1070601852    0.093048128
## #Opnn#RmFrDbt              0.008742714   0.0115740741    0.010695187
## Styls##Fshn                0.021648626   0.0086805556    0.008021390
## Bsnss#BsnssDy#SmllBsnss    0.020815987   0.0231481481    0.021925134
## myOthr                     0.006869276   0.0028935185    0.002673797
## Trvl#Trvl#                 0.017277269   0.0196759259    0.018716578
## Cltr##                              NA   0.0005787037    0.037433155
## Frgn#Wrld#AsPcfc           0.031223980   0.0306712963    0.029946524
## #Mltmd#                    0.019150708   0.0283564815    0.027807487
## TStyl##                    0.129683597   0.0584490741    0.056149733
## #U.S.#Edctn                0.050582848   0.0474537037    0.047593583
## Frgn#Wrld#                 0.026644463   0.0254629630    0.025133690
##                         err.abs.fit.sum err.abs.fit.mean .n.fit
## OpEd#Opnn#                   169.920774       0.38883472    437
## #Opnn#ThPblcEdtr               7.263079       0.45394245     16
## Styls#U.S.#                   62.317403       0.49068821    127
## Bsnss#Crsswrds/Gms#           37.637501       0.35845239    105
## Scnc#Hlth#                    67.219902       0.45418853    148
## Bsnss#Tchnlgy#                46.871671       0.22005479    213
## ##                           132.756111       0.14540647    913
## Bsnss#BsnssDy#Dlbk            96.765300       0.15383991    629
## Mtr#N.Y./Rgn#                 19.934891       0.15574134    128
## Cltr#Arts#                    60.231977       0.12292240    490
## #Opnn#RmFrDbt                  6.538023       0.15566721     42
## Styls##Fshn                    8.953809       0.08609432    104
## Bsnss#BsnssDy#SmllBsnss       13.032416       0.13032416    100
## myOthr                         3.723231       0.11282519     33
## Trvl#Trvl#                     6.800950       0.08193916     83
## Cltr##                               NA               NA     NA
## Frgn#Wrld#AsPcfc              15.573858       0.10382572    150
## #Mltmd#                        8.371241       0.09099175     92
## TStyl##                       43.550981       0.06990527    623
## #U.S.#Edctn                   15.497987       0.06377772    243
## Frgn#Wrld#                     8.953551       0.06994962    128
##                         err.abs.OOB.sum err.abs.OOB.mean
## OpEd#Opnn#                   46.5737830       0.52330093
## #Opnn#ThPblcEdtr              2.0080675       0.50201689
## Styls#U.S.#                  23.6328413       0.47265683
## Bsnss#Crsswrds/Gms#           8.3911301       0.46617390
## Scnc#Hlth#                   22.1189167       0.46081076
## Bsnss#Tchnlgy#               29.7319043       0.23596749
## ##                           78.0147094       0.21028224
## Bsnss#BsnssDy#Dlbk           66.1713471       0.20486485
## Mtr#N.Y./Rgn#                13.4746063       0.19249438
## Cltr#Arts#                   34.6944525       0.18753758
## #Opnn#RmFrDbt                 3.7425266       0.18712633
## Styls##Fshn                   2.1565235       0.14376823
## Bsnss#BsnssDy#SmllBsnss       5.6217083       0.14054271
## myOthr                        0.5675988       0.11351975
## Trvl#Trvl#                    3.6251966       0.10662343
## Cltr##                        0.1032201       0.10322006
## Frgn#Wrld#AsPcfc              5.3470704       0.10088812
## #Mltmd#                       4.8359493       0.09869284
## TStyl##                       9.5035315       0.09409437
## #U.S.#Edctn                   6.0212009       0.07342928
## Frgn#Wrld#                    3.1996976       0.07272040
##           .n.OOB           .n.Fit           .n.Tst   .freqRatio.Fit 
##      1728.000000               NA      1870.000000               NA 
##   .freqRatio.OOB   .freqRatio.Tst  err.abs.fit.sum err.abs.fit.mean 
##         1.000000         1.000000               NA               NA 
##           .n.fit  err.abs.OOB.sum err.abs.OOB.mean 
##               NA       369.535982         4.690731
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
## 1 fit.models_2_bgn          1          0    teardown 276.531  NA      NA
```

```r
glb_chunks_df <- myadd_chunk(glb_chunks_df, "fit.models", major.inc=FALSE)
```

```
##         label step_major step_minor label_minor     bgn     end elapsed
## 12 fit.models          6          2           2 261.729 276.542  14.813
## 13 fit.models          6          3           3 276.543      NA      NA
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
## 13        fit.models          6          3           3 276.543 282.563
## 14 fit.data.training          7          0           0 282.564      NA
##    elapsed
## 13    6.02
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
##          8   0.8736 0.44378   0.003353 0.02080         
##         16   0.9015 0.63767   0.006518 0.02345         
##         32   0.9010 0.63663   0.006363 0.02247         
##         55   0.9029 0.64631   0.006957 0.02433        *
## 
## The top 5 variables (out of 55):
##    WordCount.log1p, WordCount.root2, WordCount.nexp, NDSSName.my.fctrOpEd#Opnn#, PubDate.day.minutes.poly.1
## 
##  [1] "WordCount.log1p"                        
##  [2] "WordCount.root2"                        
##  [3] "WordCount.nexp"                         
##  [4] "NDSSName.my.fctrOpEd#Opnn#"             
##  [5] "PubDate.day.minutes.poly.1"             
##  [6] "PubDate.day.minutes.poly.4"             
##  [7] "PubDate.hour.fctr(15.3,23]"             
##  [8] "PubDate.last2.log1p"                    
##  [9] "PubDate.last4.log1p"                    
## [10] "NDSSName.my.fctrScnc#Hlth#"             
## [11] "NDSSName.my.fctrBsnss#Crsswrds/Gms#"    
## [12] "PubDate.day.minutes.poly.5"             
## [13] "PubDate.last8.log1p"                    
## [14] "NDSSName.my.fctrStyls#U.S.#"            
## [15] "PubDate.wkend"                          
## [16] "PubDate.last16.log1p"                   
## [17] "PubDate.day.minutes.poly.2"             
## [18] "PubDate.juliandate"                     
## [19] "PubDate.wkday.fctr6"                    
## [20] "PubDate.month.fctr11"                   
## [21] "PubDate.second.fctr(14.8,29.5]"         
## [22] "PubDate.date.fctr(7,13]"                
## [23] "PubDate.last32.log1p"                   
## [24] ".rnorm"                                 
## [25] "PubDate.wkday.fctr1"                    
## [26] "PubDate.day.minutes.poly.3"             
## [27] "PubDate.date.fctr(25,31]"               
## [28] "PubDate.hour.fctr(7.67,15.3]"           
## [29] "PubDate.minute.fctr(14.8,29.5]"         
## [30] "PubDate.month.fctr10"                   
## [31] "NDSSName.my.fctrBsnss#Tchnlgy#"         
## [32] "PubDate.second.fctr(29.5,44.2]"         
## [33] "NDSSName.my.fctrmyOthr"                 
## [34] "PubDate.date.fctr(13,19]"               
## [35] "PubDate.wkday.fctr3"                    
## [36] "PubDate.minute.fctr(44.2,59.1]"         
## [37] "PubDate.wkday.fctr4"                    
## [38] "PubDate.second.fctr(44.2,59.1]"         
## [39] "NDSSName.my.fctr#Opnn#RmFrDbt"          
## [40] "PubDate.date.fctr(19,25]"               
## [41] "NDSSName.my.fctrMtr#N.Y./Rgn#"          
## [42] "NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss"
## [43] "NDSSName.my.fctrTrvl#Trvl#"             
## [44] "NDSSName.my.fctrStyls##Fshn"            
## [45] "NDSSName.my.fctr#Mltmd#"                
## [46] "PubDate.wkday.fctr2"                    
## [47] "NDSSName.my.fctrFrgn#Wrld#"             
## [48] "NDSSName.my.fctrFrgn#Wrld#AsPcfc"       
## [49] "PubDate.wkday.fctr5"                    
## [50] "PubDate.minute.fctr(29.5,44.2]"         
## [51] "NDSSName.my.fctr#U.S.#Edctn"            
## [52] "NDSSName.my.fctrCltr#Arts#"             
## [53] "NDSSName.my.fctrBsnss#BsnssDy#Dlbk"     
## [54] "NDSSName.my.fctr##"                     
## [55] "NDSSName.my.fctrTStyl##"
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.data.training_0-1.png) 

```
## [1] "fitting model: Final##rcv#glmnet"
## [1] "    indep_vars: WordCount.root2,WordCount.log1p,NDSSName.my.fctr,PubDate.day.minutes.poly.1,PubDate.hour.fctr,PubDate.wkend,PubDate.day.minutes.poly.4,PubDate.day.minutes.poly.2,PubDate.last4.log1p,PubDate.last2.log1p,PubDate.last8.log1p,PubDate.last16.log1p,PubDate.day.minutes.poly.3,PubDate.last32.log1p,PubDate.month.fctr,PubDate.juliandate,.rnorm,PubDate.date.fctr,PubDate.second.fctr,PubDate.minute.fctr,PubDate.wkday.fctr,WordCount.nexp,PubDate.day.minutes.poly.5"
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
## Fitting alpha = 0.55, lambda = 0.00361 on full training set
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.data.training_0-2.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.data.training_0-3.png) 

```
##             Length Class      Mode     
## a0            93   -none-     numeric  
## beta        5301   dgCMatrix  S4       
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
## xNames        57   -none-     character
## problemType    1   -none-     character
## tuneValue      2   data.frame list     
## obsLevels      2   -none-     character
## [1] "min lambda > lambdaOpt:"
##                             (Intercept) 
##                             -6.71740932 
##                 NDSSName.my.fctr#Mltmd# 
##                             -0.95451508 
##           NDSSName.my.fctr#Opnn#RmFrDbt 
##                             -3.04595818 
##        NDSSName.my.fctr#Opnn#ThPblcEdtr 
##                              3.15684109 
##             NDSSName.my.fctr#U.S.#Edctn 
##                             -1.88518827 
##      NDSSName.my.fctrBsnss#BsnssDy#Dlbk 
##                             -0.24938082 
## NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss 
##                             -0.73887006 
##     NDSSName.my.fctrBsnss#Crsswrds/Gms# 
##                              3.21815213 
##          NDSSName.my.fctrBsnss#Tchnlgy# 
##                              0.44788750 
##              NDSSName.my.fctrCltr#Arts# 
##                             -0.18319632 
##              NDSSName.my.fctrFrgn#Wrld# 
##                             -1.20835894 
##        NDSSName.my.fctrFrgn#Wrld#AsPcfc 
##                             -1.61909013 
##              NDSSName.my.fctrOpEd#Opnn# 
##                              3.71825778 
##              NDSSName.my.fctrScnc#Hlth# 
##                              2.69372416 
##             NDSSName.my.fctrStyls##Fshn 
##                             -1.34904273 
##             NDSSName.my.fctrStyls#U.S.# 
##                              2.46572272 
##                 NDSSName.my.fctrTStyl## 
##                             -1.33936830 
##              NDSSName.my.fctrTrvl#Trvl# 
##                             -0.68785697 
##                  NDSSName.my.fctrmyOthr 
##                             -1.37260505 
##                 PubDate.date.fctr(7,13] 
##                              0.03756858 
##                PubDate.date.fctr(13,19] 
##                             -0.04997442 
##                PubDate.date.fctr(25,31] 
##                              0.02915688 
##              PubDate.day.minutes.poly.1 
##                             14.63604197 
##              PubDate.day.minutes.poly.2 
##                             19.28093122 
##              PubDate.day.minutes.poly.3 
##                              3.42079095 
##              PubDate.day.minutes.poly.4 
##                              2.07653730 
##            PubDate.hour.fctr(7.67,15.3] 
##                              0.24042550 
##                    PubDate.last16.log1p 
##                              0.06699389 
##                     PubDate.last2.log1p 
##                              0.01463237 
##          PubDate.minute.fctr(29.5,44.2] 
##                             -0.15908416 
##                    PubDate.month.fctr11 
##                             -0.07224476 
##          PubDate.second.fctr(44.2,59.1] 
##                             -0.08388152 
##                     PubDate.wkday.fctr1 
##                              0.09440272 
##                     PubDate.wkday.fctr2 
##                             -0.03634627 
##                     PubDate.wkday.fctr5 
##                             -0.15064565 
##                     PubDate.wkday.fctr6 
##                             -0.21128855 
##                           PubDate.wkend 
##                              0.42922277 
##                         WordCount.log1p 
##                              0.38451519 
##                         WordCount.root2 
##                              0.05236843 
## [1] "max lambda < lambdaOpt:"
##                             (Intercept) 
##                           -6.830834e+00 
##                 NDSSName.my.fctr#Mltmd# 
##                           -1.020841e+00 
##           NDSSName.my.fctr#Opnn#RmFrDbt 
##                           -3.136910e+00 
##        NDSSName.my.fctr#Opnn#ThPblcEdtr 
##                            3.164358e+00 
##             NDSSName.my.fctr#U.S.#Edctn 
##                           -1.985122e+00 
##      NDSSName.my.fctrBsnss#BsnssDy#Dlbk 
##                           -2.726760e-01 
## NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss 
##                           -7.755608e-01 
##     NDSSName.my.fctrBsnss#Crsswrds/Gms# 
##                            3.214818e+00 
##          NDSSName.my.fctrBsnss#Tchnlgy# 
##                            4.430681e-01 
##              NDSSName.my.fctrCltr#Arts# 
##                           -2.116312e-01 
##              NDSSName.my.fctrFrgn#Wrld# 
##                           -1.304342e+00 
##        NDSSName.my.fctrFrgn#Wrld#AsPcfc 
##                           -1.682973e+00 
##              NDSSName.my.fctrOpEd#Opnn# 
##                            3.721324e+00 
##              NDSSName.my.fctrScnc#Hlth# 
##                            2.688810e+00 
##             NDSSName.my.fctrStyls##Fshn 
##                           -1.425434e+00 
##             NDSSName.my.fctrStyls#U.S.# 
##                            2.459984e+00 
##                 NDSSName.my.fctrTStyl## 
##                           -1.387084e+00 
##              NDSSName.my.fctrTrvl#Trvl# 
##                           -7.480692e-01 
##                  NDSSName.my.fctrmyOthr 
##                           -1.482128e+00 
##                 PubDate.date.fctr(7,13] 
##                            4.406889e-02 
##                PubDate.date.fctr(13,19] 
##                           -5.356303e-02 
##                PubDate.date.fctr(25,31] 
##                            3.547036e-02 
##              PubDate.day.minutes.poly.1 
##                            1.496003e+01 
##              PubDate.day.minutes.poly.2 
##                            2.037019e+01 
##              PubDate.day.minutes.poly.3 
##                            3.488587e+00 
##              PubDate.day.minutes.poly.4 
##                            1.513324e+00 
##            PubDate.hour.fctr(7.67,15.3] 
##                            2.728426e-01 
##                      PubDate.juliandate 
##                           -1.676444e-05 
##                    PubDate.last16.log1p 
##                            7.368211e-02 
##                     PubDate.last2.log1p 
##                            1.445898e-02 
##          PubDate.minute.fctr(29.5,44.2] 
##                           -1.653137e-01 
##                    PubDate.month.fctr11 
##                           -7.889861e-02 
##          PubDate.second.fctr(44.2,59.1] 
##                           -9.011477e-02 
##                     PubDate.wkday.fctr1 
##                            9.895137e-02 
##                     PubDate.wkday.fctr2 
##                           -4.419254e-02 
##                     PubDate.wkday.fctr5 
##                           -1.577937e-01 
##                     PubDate.wkday.fctr6 
##                           -2.387115e-01 
##                           PubDate.wkend 
##                            4.346004e-01 
##                         WordCount.log1p 
##                            3.934762e-01 
##                         WordCount.root2 
##                            5.234272e-02
```

![](NYTBlogs3_feat_PubDate_files/figure-html/fit.data.training_0-4.png) ![](NYTBlogs3_feat_PubDate_files/figure-html/fit.data.training_0-5.png) 

```
##          Prediction
## Reference    N    Y
##         N 4967  472
##         Y  220  873
##       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull 
##   8.940600e-01   6.518908e-01   8.863433e-01   9.014229e-01   8.326699e-01 
## AccuracyPValue  McnemarPValue 
##   3.486892e-45   1.406613e-21 
##                  id
## 1 Final##rcv#glmnet
##                                                                                                                                                                                                                                                                                                                                                                                                                                                                    feats
## 1 WordCount.root2,WordCount.log1p,NDSSName.my.fctr,PubDate.day.minutes.poly.1,PubDate.hour.fctr,PubDate.wkend,PubDate.day.minutes.poly.4,PubDate.day.minutes.poly.2,PubDate.last4.log1p,PubDate.last2.log1p,PubDate.last8.log1p,PubDate.last16.log1p,PubDate.day.minutes.poly.3,PubDate.last32.log1p,PubDate.month.fctr,PubDate.juliandate,.rnorm,PubDate.date.fctr,PubDate.second.fctr,PubDate.minute.fctr,PubDate.wkday.fctr,WordCount.nexp,PubDate.day.minutes.poly.5
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1              25                     37.175                 0.898
##   max.AUCpROC.fit max.Sens.fit max.Spec.fit max.AUCROCR.fit
## 1       0.8012801    0.9612061    0.6413541       0.9361337
##   opt.prob.threshold.fit max.f.score.fit max.Accuracy.fit
## 1                    0.2       0.7161608        0.9063073
##   max.AccuracyLower.fit max.AccuracyUpper.fit max.Kappa.fit
## 1             0.8863433             0.9014229     0.6394839
##   max.AccuracySD.fit max.KappaSD.fit
## 1        0.005622392      0.02614856
```

```r
rm(ret_lst)
glb_chunks_df <- myadd_chunk(glb_chunks_df, "fit.data.training", major.inc=FALSE)
```

```
##                label step_major step_minor label_minor     bgn     end
## 14 fit.data.training          7          0           0 282.564 347.894
## 15 fit.data.training          7          1           1 347.894      NA
##    elapsed
## 14   65.33
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
##                                         All.X##rcv#glmnet.imp        imp
## PubDate.day.minutes.poly.2                          22.879214 100.000000
## PubDate.day.minutes.poly.1                         100.000000  78.090459
## NDSSName.my.fctrOpEd#Opnn#                          29.608813  29.735703
## PubDate.day.minutes.poly.3                           6.181073  28.574510
## NDSSName.my.fctrBsnss#Crsswrds/Gms#                 27.135657  27.538396
## NDSSName.my.fctr#Opnn#ThPblcEdtr                    24.849247  27.293763
## NDSSName.my.fctrScnc#Hlth#                          24.964890  25.245138
## NDSSName.my.fctrStyls#U.S.#                         23.449567  24.247822
## PubDate.day.minutes.poly.4                          46.118345  21.362667
## NDSSName.my.fctrBsnss#Tchnlgy#                       6.181073  15.438969
## PubDate.wkend                                        7.555716  15.379160
## WordCount.log1p                                      7.599061  15.191568
## PubDate.hour.fctr(7.67,15.3]                         6.181073  14.612298
## PubDate.wkday.fctr1                                  6.181073  13.915416
## PubDate.last16.log1p                                 6.181073  13.800288
## WordCount.root2                                      6.405483  13.722145
## PubDate.date.fctr(7,13]                              6.181073  13.671403
## PubDate.date.fctr(25,31]                             6.181073  13.634277
## PubDate.last2.log1p                                  6.307833  13.557058
## .rnorm                                               6.181073  13.493535
## NDSSName.my.fctrCltr##                               6.181073  13.493535
## NDSSName.my.fctrMtr#N.Y./Rgn#                        6.181073  13.493535
## PubDate.date.fctr(19,25]                             6.181073  13.493535
## PubDate.day.minutes.poly.5                           6.181073  13.493535
## PubDate.hour.fctr(15.3,23]                           6.568901  13.493535
## PubDate.last32.log1p                                 6.181073  13.493535
## PubDate.last4.log1p                                  6.403279  13.493535
## PubDate.last8.log1p                                  6.231915  13.493535
## PubDate.minute.fctr(14.8,29.5]                       6.181073  13.493535
## PubDate.minute.fctr(44.2,59.1]                       6.181073  13.493535
## PubDate.month.fctr10                                 6.181073  13.493535
## PubDate.month.fctr12                                 6.181073  13.493535
## PubDate.second.fctr(14.8,29.5]                       6.181073  13.493535
## PubDate.second.fctr(29.5,44.2]                       6.181073  13.493535
## PubDate.wkday.fctr3                                  6.181073  13.493535
## PubDate.wkday.fctr4                                  6.181073  13.493535
## WordCount.nexp                                       6.181073  13.493535
## PubDate.juliandate                                   6.181073  13.493499
## PubDate.wkday.fctr2                                  6.181073  13.318140
## PubDate.date.fctr(13,19]                             6.181073  13.267690
## PubDate.month.fctr11                                 6.181073  13.163927
## PubDate.second.fctr(44.2,59.1]                       6.181073  13.114010
## PubDate.wkday.fctr5                                  6.181073  12.820541
## PubDate.minute.fctr(29.5,44.2]                       6.181073  12.785649
## NDSSName.my.fctrCltr#Arts#                           4.494666  12.633131
## PubDate.wkday.fctr6                                  6.181073  12.512620
## NDSSName.my.fctrBsnss#BsnssDy#Dlbk                   4.575817  12.355071
## NDSSName.my.fctrTrvl#Trvl#                           4.901316  10.361955
## NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss              4.791503  10.189240
## NDSSName.my.fctr#Mltmd#                              5.672845   9.184598
## NDSSName.my.fctrFrgn#Wrld#                           4.664559   8.013114
## NDSSName.my.fctrTStyl##                              2.131448   7.543734
## NDSSName.my.fctrStyls##Fshn                          3.582132   7.440496
## NDSSName.my.fctrmyOthr                               6.019416   7.267138
## NDSSName.my.fctrFrgn#Wrld#AsPcfc                     3.104915   6.287949
## NDSSName.my.fctr#U.S.#Edctn                          3.501802   5.049357
## NDSSName.my.fctr#Opnn#RmFrDbt                        0.000000   0.000000
##                                         Final##rcv#glmnet.imp
## PubDate.day.minutes.poly.2                         100.000000
## PubDate.day.minutes.poly.1                          78.090459
## NDSSName.my.fctrOpEd#Opnn#                          29.735703
## PubDate.day.minutes.poly.3                          28.574510
## NDSSName.my.fctrBsnss#Crsswrds/Gms#                 27.538396
## NDSSName.my.fctr#Opnn#ThPblcEdtr                    27.293763
## NDSSName.my.fctrScnc#Hlth#                          25.245138
## NDSSName.my.fctrStyls#U.S.#                         24.247822
## PubDate.day.minutes.poly.4                          21.362667
## NDSSName.my.fctrBsnss#Tchnlgy#                      15.438969
## PubDate.wkend                                       15.379160
## WordCount.log1p                                     15.191568
## PubDate.hour.fctr(7.67,15.3]                        14.612298
## PubDate.wkday.fctr1                                 13.915416
## PubDate.last16.log1p                                13.800288
## WordCount.root2                                     13.722145
## PubDate.date.fctr(7,13]                             13.671403
## PubDate.date.fctr(25,31]                            13.634277
## PubDate.last2.log1p                                 13.557058
## .rnorm                                              13.493535
## NDSSName.my.fctrCltr##                              13.493535
## NDSSName.my.fctrMtr#N.Y./Rgn#                       13.493535
## PubDate.date.fctr(19,25]                            13.493535
## PubDate.day.minutes.poly.5                          13.493535
## PubDate.hour.fctr(15.3,23]                          13.493535
## PubDate.last32.log1p                                13.493535
## PubDate.last4.log1p                                 13.493535
## PubDate.last8.log1p                                 13.493535
## PubDate.minute.fctr(14.8,29.5]                      13.493535
## PubDate.minute.fctr(44.2,59.1]                      13.493535
## PubDate.month.fctr10                                13.493535
## PubDate.month.fctr12                                13.493535
## PubDate.second.fctr(14.8,29.5]                      13.493535
## PubDate.second.fctr(29.5,44.2]                      13.493535
## PubDate.wkday.fctr3                                 13.493535
## PubDate.wkday.fctr4                                 13.493535
## WordCount.nexp                                      13.493535
## PubDate.juliandate                                  13.493499
## PubDate.wkday.fctr2                                 13.318140
## PubDate.date.fctr(13,19]                            13.267690
## PubDate.month.fctr11                                13.163927
## PubDate.second.fctr(44.2,59.1]                      13.114010
## PubDate.wkday.fctr5                                 12.820541
## PubDate.minute.fctr(29.5,44.2]                      12.785649
## NDSSName.my.fctrCltr#Arts#                          12.633131
## PubDate.wkday.fctr6                                 12.512620
## NDSSName.my.fctrBsnss#BsnssDy#Dlbk                  12.355071
## NDSSName.my.fctrTrvl#Trvl#                          10.361955
## NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss             10.189240
## NDSSName.my.fctr#Mltmd#                              9.184598
## NDSSName.my.fctrFrgn#Wrld#                           8.013114
## NDSSName.my.fctrTStyl##                              7.543734
## NDSSName.my.fctrStyls##Fshn                          7.440496
## NDSSName.my.fctrmyOthr                               7.267138
## NDSSName.my.fctrFrgn#Wrld#AsPcfc                     6.287949
## NDSSName.my.fctr#U.S.#Edctn                          5.049357
## NDSSName.my.fctr#Opnn#RmFrDbt                        0.000000
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
## 2     4168         N                        0.0411419
## 3     5647         N                        0.1290645
## 4      302         N                               NA
##   Pplr.fctr.All.X..rcv.glmnet Pplr.fctr.All.X..rcv.glmnet.err
## 1                        <NA>                              NA
## 2                           N                           FALSE
## 3                           Y                            TRUE
## 4                        <NA>                              NA
##   Pplr.fctr.All.X..rcv.glmnet.err.abs Pplr.fctr.All.X..rcv.glmnet.is.acc
## 1                                  NA                                 NA
## 2                           0.0411419                               TRUE
## 3                           0.1290645                              FALSE
## 4                                  NA                                 NA
##   Pplr.fctr.Final..rcv.glmnet.prob Pplr.fctr.Final..rcv.glmnet
## 1                      0.022283474                           N
## 2                      0.007666321                           N
## 3                      0.105843562                           Y
## 4                      0.436516680                           Y
##   Pplr.fctr.Final..rcv.glmnet.err Pplr.fctr.Final..rcv.glmnet.err.abs
## 1                           FALSE                         0.022283474
## 2                           FALSE                         0.007666321
## 3                            TRUE                         0.105843562
## 4                            TRUE                         0.436516680
##   Pplr.fctr.Final..rcv.glmnet.is.acc Pplr.fctr.Final..rcv.glmnet.accurate
## 1                               TRUE                                 TRUE
## 2                               TRUE                                 TRUE
## 3                              FALSE                                FALSE
## 4                              FALSE                                FALSE
##   Pplr.fctr.Final..rcv.glmnet.error .label
## 1                       0.000000000   1065
## 2                       0.000000000   4168
## 3                       0.005843562   5647
## 4                       0.336516680    302
## [1] "Inaccurate: "
##   UniqueID Pplr.fctr Pplr.fctr.All.X..rcv.glmnet.prob
## 1     2182         Y                       0.04423435
## 2     4352         Y                               NA
## 3     5486         Y                               NA
## 4     4721         Y                       0.06779771
## 5     1696         Y                       0.05294423
## 6      364         Y                       0.06987153
##   Pplr.fctr.All.X..rcv.glmnet Pplr.fctr.All.X..rcv.glmnet.err
## 1                           N                            TRUE
## 2                        <NA>                              NA
## 3                        <NA>                              NA
## 4                           N                            TRUE
## 5                           N                            TRUE
## 6                           N                            TRUE
##   Pplr.fctr.All.X..rcv.glmnet.err.abs Pplr.fctr.All.X..rcv.glmnet.is.acc
## 1                           0.9557656                              FALSE
## 2                                  NA                                 NA
## 3                                  NA                                 NA
## 4                           0.9322023                              FALSE
## 5                           0.9470558                              FALSE
## 6                           0.9301285                              FALSE
##   Pplr.fctr.Final..rcv.glmnet.prob Pplr.fctr.Final..rcv.glmnet
## 1                      0.006646728                           N
## 2                      0.009472057                           N
## 3                      0.013199051                           N
## 4                      0.013513861                           N
## 5                      0.015668226                           N
## 6                      0.019168273                           N
##   Pplr.fctr.Final..rcv.glmnet.err Pplr.fctr.Final..rcv.glmnet.err.abs
## 1                            TRUE                           0.9933533
## 2                            TRUE                           0.9905279
## 3                            TRUE                           0.9868009
## 4                            TRUE                           0.9864861
## 5                            TRUE                           0.9843318
## 6                            TRUE                           0.9808317
##   Pplr.fctr.Final..rcv.glmnet.is.acc Pplr.fctr.Final..rcv.glmnet.accurate
## 1                              FALSE                                FALSE
## 2                              FALSE                                FALSE
## 3                              FALSE                                FALSE
## 4                              FALSE                                FALSE
## 5                              FALSE                                FALSE
## 6                              FALSE                                FALSE
##   Pplr.fctr.Final..rcv.glmnet.error
## 1                       -0.09335327
## 2                       -0.09052794
## 3                       -0.08680095
## 4                       -0.08648614
## 5                       -0.08433177
## 6                       -0.08083173
##      UniqueID Pplr.fctr Pplr.fctr.All.X..rcv.glmnet.prob
## 190       789         N                        0.1059522
## 220      3588         N                               NA
## 380      3380         N                               NA
## 410        72         N                        0.1364144
## 905      2226         N                        0.1913497
## 1105     5462         N                               NA
##      Pplr.fctr.All.X..rcv.glmnet Pplr.fctr.All.X..rcv.glmnet.err
## 190                            Y                            TRUE
## 220                         <NA>                              NA
## 380                         <NA>                              NA
## 410                            Y                            TRUE
## 905                            Y                            TRUE
## 1105                        <NA>                              NA
##      Pplr.fctr.All.X..rcv.glmnet.err.abs
## 190                            0.1059522
## 220                                   NA
## 380                                   NA
## 410                            0.1364144
## 905                            0.1913497
## 1105                                  NA
##      Pplr.fctr.All.X..rcv.glmnet.is.acc Pplr.fctr.Final..rcv.glmnet.prob
## 190                               FALSE                        0.1129611
## 220                                  NA                        0.1179667
## 380                                  NA                        0.1378456
## 410                               FALSE                        0.1413828
## 905                               FALSE                        0.3794899
## 1105                                 NA                        0.7300482
##      Pplr.fctr.Final..rcv.glmnet Pplr.fctr.Final..rcv.glmnet.err
## 190                            Y                            TRUE
## 220                            Y                            TRUE
## 380                            Y                            TRUE
## 410                            Y                            TRUE
## 905                            Y                            TRUE
## 1105                           Y                            TRUE
##      Pplr.fctr.Final..rcv.glmnet.err.abs
## 190                            0.1129611
## 220                            0.1179667
## 380                            0.1378456
## 410                            0.1413828
## 905                            0.3794899
## 1105                           0.7300482
##      Pplr.fctr.Final..rcv.glmnet.is.acc
## 190                               FALSE
## 220                               FALSE
## 380                               FALSE
## 410                               FALSE
## 905                               FALSE
## 1105                              FALSE
##      Pplr.fctr.Final..rcv.glmnet.accurate
## 190                                 FALSE
## 220                                 FALSE
## 380                                 FALSE
## 410                                 FALSE
## 905                                 FALSE
## 1105                                FALSE
##      Pplr.fctr.Final..rcv.glmnet.error
## 190                         0.01296111
## 220                         0.01796670
## 380                         0.03784555
## 410                         0.04138278
## 905                         0.27948986
## 1105                        0.63004823
##      UniqueID Pplr.fctr Pplr.fctr.All.X..rcv.glmnet.prob
## 1182      770         N                               NA
## 1183     2179         N                        0.7720249
## 1184      472         N                               NA
## 1185     2995         N                               NA
## 1186     1612         N                        0.7937303
## 1187     1448         N                               NA
##      Pplr.fctr.All.X..rcv.glmnet Pplr.fctr.All.X..rcv.glmnet.err
## 1182                        <NA>                              NA
## 1183                           Y                            TRUE
## 1184                        <NA>                              NA
## 1185                        <NA>                              NA
## 1186                           Y                            TRUE
## 1187                        <NA>                              NA
##      Pplr.fctr.All.X..rcv.glmnet.err.abs
## 1182                                  NA
## 1183                           0.7720249
## 1184                                  NA
## 1185                                  NA
## 1186                           0.7937303
## 1187                                  NA
##      Pplr.fctr.All.X..rcv.glmnet.is.acc Pplr.fctr.Final..rcv.glmnet.prob
## 1182                                 NA                        0.9518142
## 1183                              FALSE                        0.9542775
## 1184                                 NA                        0.9580719
## 1185                                 NA                        0.9628305
## 1186                              FALSE                        0.9643315
## 1187                                 NA                        0.9689304
##      Pplr.fctr.Final..rcv.glmnet Pplr.fctr.Final..rcv.glmnet.err
## 1182                           Y                            TRUE
## 1183                           Y                            TRUE
## 1184                           Y                            TRUE
## 1185                           Y                            TRUE
## 1186                           Y                            TRUE
## 1187                           Y                            TRUE
##      Pplr.fctr.Final..rcv.glmnet.err.abs
## 1182                           0.9518142
## 1183                           0.9542775
## 1184                           0.9580719
## 1185                           0.9628305
## 1186                           0.9643315
## 1187                           0.9689304
##      Pplr.fctr.Final..rcv.glmnet.is.acc
## 1182                              FALSE
## 1183                              FALSE
## 1184                              FALSE
## 1185                              FALSE
## 1186                              FALSE
## 1187                              FALSE
##      Pplr.fctr.Final..rcv.glmnet.accurate
## 1182                                FALSE
## 1183                                FALSE
## 1184                                FALSE
## 1185                                FALSE
## 1186                                FALSE
## 1187                                FALSE
##      Pplr.fctr.Final..rcv.glmnet.error
## 1182                         0.8518142
## 1183                         0.8542775
## 1184                         0.8580719
## 1185                         0.8628305
## 1186                         0.8643315
## 1187                         0.8689304
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
## 15 fit.data.training          7          1           1 347.894 358.247
## 16  predict.data.new          8          0           0 358.247      NA
##    elapsed
## 15  10.353
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
##    NDSSName.my.fctr .n.Trn.N .n.Trn.Y .n.New.N .n.New.Y
## 5       #U.S.#Edctn      325       NA       87        2
## 10           Cltr##        1       NA       47       23
## .n.Trn.N .n.Trn.Y .n.New.N .n.New.Y 
##      326        0      134       25
```

```
## Loading required package: tidyr
## 
## Attaching package: 'tidyr'
## 
## The following object is masked from 'package:Matrix':
## 
##     expand
## 
## The following object is masked from 'package:mice':
## 
##     complete
```

```
## [1] "OOBobs Pplr.fctr.All.X..rcv.glmnet Y: min < min of Train range: 12"
##      UniqueID Pplr.fctr.All.X..rcv.glmnet PubDate.day.minutes.poly.2
## 1751     1751                           Y               -0.001148534
## 639       639                           Y                0.009120600
## 1542     1542                           Y               -0.008640045
## 1535     1535                           Y               -0.008325008
## 2645     2645                           Y                0.012594250
## 697       697                           Y               -0.008337839
## 2382     2382                           Y               -0.006921772
## 5630     5630                           Y               -0.007831741
## 1677     1677                           Y               -0.008708358
## 6435     6435                           Y                0.002105741
## 1871     1871                           Y               -0.001039862
## 4223     4223                           Y               -0.008758791
##      PubDate.day.minutes.poly.4 PubDate.last2.log1p PubDate.last32.log1p
## 1751               -0.007894421            5.463832             9.151227
## 639                -0.017000978            2.995732            10.597135
## 1542                0.008801265            6.313548             9.217117
## 1535                0.007631110            5.111988             9.222862
## 2645               -0.018240126            8.943637            10.695801
## 697                 0.007674714            6.338594             9.641213
## 2382                0.003523823            6.318968             9.397235
## 5630                0.006070276            6.104793             9.495294
## 1677                0.009654246            4.890349             9.758866
## 6435               -0.010140008            7.624131            10.866967
## 1871               -0.008047440            4.804021             9.176990
## 4223                0.009518864            8.350902            11.242625
##      PubDate.last8.log1p WordCount.log1p WordCount.root2
## 1751            7.571988        7.070724       34.292856
## 639            10.093116        6.986566       32.878564
## 1542            7.846981        7.135687       35.425979
## 1535            7.569928        6.463029       25.298221
## 2645           10.083765        6.635947       27.586228
## 697             7.019297        6.415097       24.698178
## 2382            7.138867        5.755742       17.748239
## 5630            7.152269        3.737670        6.403124
## 1677            7.124478        6.756932       29.308702
## 6435            8.695172        0.000000        0.000000
## 1871            7.187657        7.034388       33.674916
## 4223           10.840972        7.274480       37.973675
##                                                    id      cor.y
## PubDate.day.minutes.poly.2 PubDate.day.minutes.poly.2 0.07097772
## PubDate.day.minutes.poly.4 PubDate.day.minutes.poly.4 0.07394139
## PubDate.last2.log1p               PubDate.last2.log1p 0.06567954
## PubDate.last32.log1p             PubDate.last32.log1p 0.02254251
## PubDate.last8.log1p               PubDate.last8.log1p 0.05657458
## WordCount.log1p                       WordCount.log1p 0.25431963
## WordCount.root2                       WordCount.root2 0.29212068
##                            exclude.as.feat  cor.y.abs          cor.high.X
## PubDate.day.minutes.poly.2           FALSE 0.07097772                <NA>
## PubDate.day.minutes.poly.4           FALSE 0.07394139                <NA>
## PubDate.last2.log1p                  FALSE 0.06567954 PubDate.last4.log1p
## PubDate.last32.log1p                 FALSE 0.02254251                <NA>
## PubDate.last8.log1p                  FALSE 0.05657458 PubDate.last4.log1p
## WordCount.log1p                      FALSE 0.25431963     WordCount.root2
## WordCount.root2                      FALSE 0.29212068                <NA>
##                            freqRatio percentUnique zeroVar   nzv
## PubDate.day.minutes.poly.2  1.225490      18.08022   FALSE FALSE
## PubDate.day.minutes.poly.4  1.225490      18.08022   FALSE FALSE
## PubDate.last2.log1p         1.375000      51.16350   FALSE FALSE
## PubDate.last32.log1p        1.000000      91.13595   FALSE FALSE
## PubDate.last8.log1p         1.166667      75.15309   FALSE FALSE
## WordCount.log1p             2.315789      24.15799   FALSE FALSE
## WordCount.root2             2.315789      24.15799   FALSE FALSE
##                            is.cor.y.abs.low interaction.feat
## PubDate.day.minutes.poly.2            FALSE               NA
## PubDate.day.minutes.poly.4            FALSE               NA
## PubDate.last2.log1p                   FALSE               NA
## PubDate.last32.log1p                  FALSE               NA
## PubDate.last8.log1p                   FALSE               NA
## WordCount.log1p                       FALSE               NA
## WordCount.root2                       FALSE               NA
##                            shapiro.test.p.value rsp_var_raw id_var rsp_var
## PubDate.day.minutes.poly.2         8.020999e-64       FALSE     NA      NA
## PubDate.day.minutes.poly.4         1.523136e-47       FALSE     NA      NA
## PubDate.last2.log1p                4.200850e-29       FALSE     NA      NA
## PubDate.last32.log1p               1.048418e-44       FALSE     NA      NA
## PubDate.last8.log1p                1.021260e-42       FALSE     NA      NA
## WordCount.log1p                    1.576866e-49       FALSE     NA      NA
## WordCount.root2                    4.556481e-30       FALSE     NA      NA
##                                     max          min max.Pplr.fctr.N
## PubDate.day.minutes.poly.2   0.04268445 -0.008758791      0.04268445
## PubDate.day.minutes.poly.4   0.06677441 -0.018327397      0.06543120
## PubDate.last2.log1p         10.91197453  0.693147181     10.91197453
## PubDate.last32.log1p        12.32340669  8.835792367     12.30546086
## PubDate.last8.log1p         11.62246125  6.666956792     11.43577441
## WordCount.log1p              9.29771002  0.000000000      8.81966535
## WordCount.root2            104.46051886  0.000000000     82.24962006
##                            max.Pplr.fctr.Y min.Pplr.fctr.N min.Pplr.fctr.Y
## PubDate.day.minutes.poly.2      0.04254377    -0.008758791    -0.008758717
## PubDate.day.minutes.poly.4      0.06149053    -0.018327397    -0.018219595
## PubDate.last2.log1p            10.89191288     0.693147181     3.135494216
## PubDate.last32.log1p           12.17840850     9.125653564     9.237663668
## PubDate.last8.log1p            11.39428831     6.666956792     7.199678346
## WordCount.log1p                 9.29771002     0.000000000     1.945910149
## WordCount.root2               104.46051886     0.000000000     2.449489743
##                            max.Pplr.fctr.All.X..rcv.glmnet.N
## PubDate.day.minutes.poly.2                        0.04268445
## PubDate.day.minutes.poly.4                        0.06213811
## PubDate.last2.log1p                              10.71916205
## PubDate.last32.log1p                             12.17383350
## PubDate.last8.log1p                              11.40150216
## WordCount.log1p                                   7.05961763
## WordCount.root2                                  34.10278581
##                            max.Pplr.fctr.All.X..rcv.glmnet.Y
## PubDate.day.minutes.poly.2                        0.04254377
## PubDate.day.minutes.poly.4                        0.06610094
## PubDate.last2.log1p                              10.84161805
## PubDate.last32.log1p                             12.21791228
## PubDate.last8.log1p                              11.62246125
## WordCount.log1p                                   9.14088311
## WordCount.root2                                  96.58157174
##                            min.Pplr.fctr.All.X..rcv.glmnet.N
## PubDate.day.minutes.poly.2                      -0.008758672
## PubDate.day.minutes.poly.4                      -0.018326850
## PubDate.last2.log1p                              1.386294361
## PubDate.last32.log1p                             9.194210990
## PubDate.last8.log1p                              6.934397210
## WordCount.log1p                                  0.000000000
## WordCount.root2                                  0.000000000
##                            min.Pplr.fctr.All.X..rcv.glmnet.Y
## PubDate.day.minutes.poly.2                      -0.008758791
## PubDate.day.minutes.poly.4                      -0.018240126
## PubDate.last2.log1p                              2.995732274
## PubDate.last32.log1p                             9.151227107
## PubDate.last8.log1p                              7.019296654
## WordCount.log1p                                  0.000000000
## WordCount.root2                                  0.000000000
##                            max.Pplr.fctr.Final..rcv.glmnet.N
## PubDate.day.minutes.poly.2                        0.04268445
## PubDate.day.minutes.poly.4                        0.06543120
## PubDate.last2.log1p                              10.74567954
## PubDate.last32.log1p                             12.32340669
## PubDate.last8.log1p                              11.27955479
## WordCount.log1p                                   7.97384438
## WordCount.root2                                  53.87949517
##                            max.Pplr.fctr.Final..rcv.glmnet.Y
## PubDate.day.minutes.poly.2                        0.04254377
## PubDate.day.minutes.poly.4                        0.06677441
## PubDate.last2.log1p                              10.75436471
## PubDate.last32.log1p                             12.30546086
## PubDate.last8.log1p                              11.33227851
## WordCount.log1p                                   8.69232228
## WordCount.root2                                  77.17512553
##                            min.Pplr.fctr.Final..rcv.glmnet.N
## PubDate.day.minutes.poly.2                      -0.008758791
## PubDate.day.minutes.poly.4                      -0.018322678
## PubDate.last2.log1p                              2.197224577
## PubDate.last32.log1p                             8.862908295
## PubDate.last8.log1p                              7.061334367
## WordCount.log1p                                  0.000000000
## WordCount.root2                                  0.000000000
##                            min.Pplr.fctr.Final..rcv.glmnet.Y
## PubDate.day.minutes.poly.2                      -0.008758672
## PubDate.day.minutes.poly.4                      -0.018203392
## PubDate.last2.log1p                              3.610917913
## PubDate.last32.log1p                             8.835792367
## PubDate.last8.log1p                              6.891625897
## WordCount.log1p                                  1.609437912
## WordCount.root2                                  2.000000000
## [1] "OOBobs Pplr.fctr.All.X..rcv.glmnet Y: max > max of Train range: 7"
##      UniqueID Pplr.fctr.All.X..rcv.glmnet PubDate.day.minutes.poly.1
## 5233     5233                           Y               -0.016600865
## 6528     6528                           Y                0.001809613
## 302       302                           Y                0.024722851
## 6435     6435                           Y               -0.013151170
## 6517     6517                           Y                0.020474279
## 3205     3205                           Y                0.002898990
## 6521     6521                           Y                0.013901702
##      PubDate.day.minutes.poly.3 PubDate.day.minutes.poly.4
## 5233                0.007145007               -0.017000978
## 6528               -0.004472838                0.007219141
## 302                 0.051829879                0.066100937
## 6435                0.010843928               -0.010140008
## 6517                0.020825101                0.010136126
## 3205               -0.005881110                0.005605493
## 6521               -0.004412350               -0.013803062
##      PubDate.day.minutes.poly.5 PubDate.last16.log1p PubDate.last32.log1p
## 5233                0.013803383             11.76273             11.81975
## 6528                0.006658270             11.95698             12.10953
## 302                 0.083442278             10.07706             10.33290
## 6435               -0.001891593             10.33864             10.86697
## 6517               -0.004924177             11.44856             12.21791
## 3205                0.008323368             11.85059             12.00544
## 6521               -0.012191759             11.68539             12.19622
##      PubDate.last8.log1p WordCount.nexp
## 5233           11.622461  5.482209e-194
## 6528           11.425547   0.000000e+00
## 302             9.741557   0.000000e+00
## 6435            8.695172   1.000000e+00
## 6517            9.753188  2.750325e-314
## 3205           11.443361   1.026188e-10
## 6521           10.414633   0.000000e+00
##                                                    id       cor.y
## PubDate.day.minutes.poly.1 PubDate.day.minutes.poly.1  0.15675348
## PubDate.day.minutes.poly.3 PubDate.day.minutes.poly.3  0.02798355
## PubDate.day.minutes.poly.4 PubDate.day.minutes.poly.4  0.07394139
## PubDate.day.minutes.poly.5 PubDate.day.minutes.poly.5 -0.05592923
## PubDate.last16.log1p             PubDate.last16.log1p  0.03845681
## PubDate.last32.log1p             PubDate.last32.log1p  0.02254251
## PubDate.last8.log1p               PubDate.last8.log1p  0.05657458
## WordCount.nexp                         WordCount.nexp -0.05320840
##                            exclude.as.feat  cor.y.abs          cor.high.X
## PubDate.day.minutes.poly.1           FALSE 0.15675348                <NA>
## PubDate.day.minutes.poly.3           FALSE 0.02798355                <NA>
## PubDate.day.minutes.poly.4           FALSE 0.07394139                <NA>
## PubDate.day.minutes.poly.5           FALSE 0.05592923                <NA>
## PubDate.last16.log1p                 FALSE 0.03845681 PubDate.last8.log1p
## PubDate.last32.log1p                 FALSE 0.02254251                <NA>
## PubDate.last8.log1p                  FALSE 0.05657458 PubDate.last4.log1p
## WordCount.nexp                       FALSE 0.05320840                <NA>
##                            freqRatio percentUnique zeroVar   nzv
## PubDate.day.minutes.poly.1  1.225490      18.08022   FALSE FALSE
## PubDate.day.minutes.poly.3  1.225490      18.08022   FALSE FALSE
## PubDate.day.minutes.poly.4  1.225490      18.08022   FALSE FALSE
## PubDate.day.minutes.poly.5  1.225490      18.08022   FALSE FALSE
## PubDate.last16.log1p        1.000000      84.50704   FALSE FALSE
## PubDate.last32.log1p        1.000000      91.13595   FALSE FALSE
## PubDate.last8.log1p         1.166667      75.15309   FALSE FALSE
## WordCount.nexp             17.761364      11.32884   FALSE FALSE
##                            is.cor.y.abs.low interaction.feat
## PubDate.day.minutes.poly.1            FALSE               NA
## PubDate.day.minutes.poly.3            FALSE               NA
## PubDate.day.minutes.poly.4            FALSE               NA
## PubDate.day.minutes.poly.5            FALSE               NA
## PubDate.last16.log1p                  FALSE               NA
## PubDate.last32.log1p                  FALSE               NA
## PubDate.last8.log1p                   FALSE               NA
## WordCount.nexp                        FALSE               NA
##                            shapiro.test.p.value rsp_var_raw id_var rsp_var
## PubDate.day.minutes.poly.1         1.590362e-18       FALSE     NA      NA
## PubDate.day.minutes.poly.3         9.822405e-52       FALSE     NA      NA
## PubDate.day.minutes.poly.4         1.523136e-47       FALSE     NA      NA
## PubDate.day.minutes.poly.5         1.157500e-41       FALSE     NA      NA
## PubDate.last16.log1p               8.465226e-47       FALSE     NA      NA
## PubDate.last32.log1p               1.048418e-44       FALSE     NA      NA
## PubDate.last8.log1p                1.021260e-42       FALSE     NA      NA
## WordCount.nexp                     9.108805e-94       FALSE     NA      NA
##                                    max         min max.Pplr.fctr.N
## PubDate.day.minutes.poly.1  0.02475916 -0.02749464      0.02468654
## PubDate.day.minutes.poly.3  0.05215301 -0.04512497      0.05150779
## PubDate.day.minutes.poly.4  0.06677441 -0.01832740      0.06543120
## PubDate.day.minutes.poly.5  0.08471756 -0.02450918      0.08217780
## PubDate.last16.log1p       11.95698288  8.00068478     11.94531808
## PubDate.last32.log1p       12.32340669  8.83579237     12.30546086
## PubDate.last8.log1p        11.62246125  6.66695679     11.43577441
## WordCount.nexp              1.00000000  0.00000000      1.00000000
##                            max.Pplr.fctr.Y min.Pplr.fctr.N min.Pplr.fctr.Y
## PubDate.day.minutes.poly.1     0.024468663     -0.02749464     -0.02745833
## PubDate.day.minutes.poly.3     0.049597025     -0.04512497     -0.04482024
## PubDate.day.minutes.poly.4     0.061490534     -0.01832740     -0.01821959
## PubDate.day.minutes.poly.5     0.074814724     -0.02450918     -0.02362780
## PubDate.last16.log1p          11.877603300      8.05452261      8.00068478
## PubDate.last32.log1p          12.178408497      9.12565356      9.23766367
## PubDate.last8.log1p           11.394288315      6.66695679      7.19967835
## WordCount.nexp                 0.002478752      0.00000000      0.00000000
##                            max.Pplr.fctr.All.X..rcv.glmnet.N
## PubDate.day.minutes.poly.1                        0.02450498
## PubDate.day.minutes.poly.3                        0.04991290
## PubDate.day.minutes.poly.4                        0.06213811
## PubDate.day.minutes.poly.5                        0.07601554
## PubDate.last16.log1p                             11.84854019
## PubDate.last32.log1p                             12.17383350
## PubDate.last8.log1p                              11.40150216
## WordCount.nexp                                    1.00000000
##                            max.Pplr.fctr.All.X..rcv.glmnet.Y
## PubDate.day.minutes.poly.1                        0.02472285
## PubDate.day.minutes.poly.3                        0.05182988
## PubDate.day.minutes.poly.4                        0.06610094
## PubDate.day.minutes.poly.5                        0.08344228
## PubDate.last16.log1p                             11.95698288
## PubDate.last32.log1p                             12.21791228
## PubDate.last8.log1p                              11.62246125
## WordCount.nexp                                    1.00000000
##                            min.Pplr.fctr.All.X..rcv.glmnet.N
## PubDate.day.minutes.poly.1                       -0.02749464
## PubDate.day.minutes.poly.3                       -0.04512497
## PubDate.day.minutes.poly.4                       -0.01832685
## PubDate.day.minutes.poly.5                       -0.02450918
## PubDate.last16.log1p                              8.23164218
## PubDate.last32.log1p                              9.19421099
## PubDate.last8.log1p                               6.93439721
## WordCount.nexp                                    0.00000000
##                            min.Pplr.fctr.All.X..rcv.glmnet.Y
## PubDate.day.minutes.poly.1                       -0.02745833
## PubDate.day.minutes.poly.3                       -0.04482024
## PubDate.day.minutes.poly.4                       -0.01824013
## PubDate.day.minutes.poly.5                       -0.02362780
## PubDate.last16.log1p                              8.29953457
## PubDate.last32.log1p                              9.15122711
## PubDate.last8.log1p                               7.01929665
## WordCount.nexp                                    0.00000000
##                            max.Pplr.fctr.Final..rcv.glmnet.N
## PubDate.day.minutes.poly.1                        0.02468654
## PubDate.day.minutes.poly.3                        0.05150779
## PubDate.day.minutes.poly.4                        0.06543120
## PubDate.day.minutes.poly.5                        0.08217780
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
## [1] "OOBobs Pplr.fctr.All.X..rcv.glmnet N: min < min of Train range: 1"
##      UniqueID Pplr.fctr.All.X..rcv.glmnet PubDate.last4.log1p
## 2383     2383                           N            4.158883
##                                      id     cor.y exclude.as.feat
## PubDate.last4.log1p PubDate.last4.log1p 0.0697764           FALSE
##                     cor.y.abs cor.high.X freqRatio percentUnique zeroVar
## PubDate.last4.log1p 0.0697764       <NA>     1.125      64.98775   FALSE
##                       nzv is.cor.y.abs.low interaction.feat
## PubDate.last4.log1p FALSE            FALSE               NA
##                     shapiro.test.p.value rsp_var_raw id_var rsp_var
## PubDate.last4.log1p         9.671475e-32       FALSE     NA      NA
##                          max      min max.Pplr.fctr.N max.Pplr.fctr.Y
## PubDate.last4.log1p 11.24276 4.158883        11.24276         11.2174
##                     min.Pplr.fctr.N min.Pplr.fctr.Y
## PubDate.last4.log1p        4.382027        5.153292
##                     max.Pplr.fctr.All.X..rcv.glmnet.N
## PubDate.last4.log1p                           10.8627
##                     max.Pplr.fctr.All.X..rcv.glmnet.Y
## PubDate.last4.log1p                          11.02225
##                     min.Pplr.fctr.All.X..rcv.glmnet.N
## PubDate.last4.log1p                          4.158883
##                     min.Pplr.fctr.All.X..rcv.glmnet.Y
## PubDate.last4.log1p                          5.897154
##                     max.Pplr.fctr.Final..rcv.glmnet.N
## PubDate.last4.log1p                          10.98656
##                     max.Pplr.fctr.Final..rcv.glmnet.Y
## PubDate.last4.log1p                          11.08721
##                     min.Pplr.fctr.Final..rcv.glmnet.N
## PubDate.last4.log1p                          4.941642
##                     min.Pplr.fctr.Final..rcv.glmnet.Y
## PubDate.last4.log1p                          5.442418
## [1] "OOBobs total range outliers: 19"
## [1] "newobs Pplr.fctr.Final..rcv.glmnet N: min < min of Train range: 29"
##      UniqueID Pplr.fctr.Final..rcv.glmnet PubDate.last32.log1p
## 7744     7744                           N             9.100860
## 7745     7745                           N             9.061608
## 7747     7747                           N             9.050172
## 7749     7749                           N             8.980550
## 7750     7750                           N             8.997518
## 7751     7751                           N             8.982310
##      UniqueID Pplr.fctr.Final..rcv.glmnet PubDate.last32.log1p
## 7745     7745                           N             9.061608
## 7761     7761                           N             9.001100
## 7763     7763                           N             9.065315
## 7766     7766                           N             9.022081
## 7936     7936                           N             8.932741
## 7941     7941                           N             8.956351
##      UniqueID Pplr.fctr.Final..rcv.glmnet PubDate.last32.log1p
## 7941     7941                           N             8.956351
## 7942     7942                           N             8.921991
## 7943     7943                           N             8.906393
## 7944     7944                           N             8.991064
## 7945     7945                           N             9.004054
## 7946     7946                           N             9.121509
##                                        id      cor.y exclude.as.feat
## PubDate.last32.log1p PubDate.last32.log1p 0.02254251           FALSE
##                       cor.y.abs cor.high.X freqRatio percentUnique zeroVar
## PubDate.last32.log1p 0.02254251       <NA>         1      91.13595   FALSE
##                        nzv is.cor.y.abs.low interaction.feat
## PubDate.last32.log1p FALSE            FALSE               NA
##                      shapiro.test.p.value rsp_var_raw id_var rsp_var
## PubDate.last32.log1p         1.048418e-44       FALSE     NA      NA
##                           max      min max.Pplr.fctr.N max.Pplr.fctr.Y
## PubDate.last32.log1p 12.32341 8.835792        12.30546        12.21791
##                      min.Pplr.fctr.N min.Pplr.fctr.Y
## PubDate.last32.log1p        9.125654         9.17699
##                      max.Pplr.fctr.All.X..rcv.glmnet.N
## PubDate.last32.log1p                          12.17383
##                      max.Pplr.fctr.All.X..rcv.glmnet.Y
## PubDate.last32.log1p                          12.21791
##                      min.Pplr.fctr.All.X..rcv.glmnet.N
## PubDate.last32.log1p                          9.194211
##                      min.Pplr.fctr.All.X..rcv.glmnet.Y
## PubDate.last32.log1p                          9.151227
##                      max.Pplr.fctr.Final..rcv.glmnet.N
## PubDate.last32.log1p                          12.32341
##                      max.Pplr.fctr.Final..rcv.glmnet.Y
## PubDate.last32.log1p                          12.30546
##                      min.Pplr.fctr.Final..rcv.glmnet.N
## PubDate.last32.log1p                          8.862908
##                      min.Pplr.fctr.Final..rcv.glmnet.Y
## PubDate.last32.log1p                          8.835792
## [1] "newobs Pplr.fctr.Final..rcv.glmnet N: max > max of Train range: 1189"
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
## 6717     6717                           N                336
## 6880     6880                           N                338
## 7332     7332                           N                345
## 7537     7537                           N                349
## 7724     7724                           N                351
## 8062     8062                           N                356
##      PubDate.last32.log1p
## 6717            10.684760
## 6880            10.930371
## 7332             9.601842
## 7537            10.000796
## 7724             9.325097
## 8062            10.218590
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
##                                        id      cor.y exclude.as.feat
## PubDate.juliandate     PubDate.juliandate 0.01436107           FALSE
## PubDate.last32.log1p PubDate.last32.log1p 0.02254251           FALSE
##                       cor.y.abs         cor.high.X freqRatio percentUnique
## PubDate.juliandate   0.01436107 PubDate.month.fctr   1.03252      1.393141
## PubDate.last32.log1p 0.02254251               <NA>   1.00000     91.135946
##                      zeroVar   nzv is.cor.y.abs.low interaction.feat
## PubDate.juliandate     FALSE FALSE            FALSE               NA
## PubDate.last32.log1p   FALSE FALSE            FALSE               NA
##                      shapiro.test.p.value rsp_var_raw id_var rsp_var
## PubDate.juliandate           1.389406e-35       FALSE     NA      NA
## PubDate.last32.log1p         1.048418e-44       FALSE     NA      NA
##                            max        min max.Pplr.fctr.N max.Pplr.fctr.Y
## PubDate.juliandate   365.00000 244.000000       334.00000       334.00000
## PubDate.last32.log1p  12.32341   8.835792        12.30546        12.21791
##                      min.Pplr.fctr.N min.Pplr.fctr.Y
## PubDate.juliandate        244.000000       244.00000
## PubDate.last32.log1p        9.125654         9.17699
##                      max.Pplr.fctr.All.X..rcv.glmnet.N
## PubDate.juliandate                           332.00000
## PubDate.last32.log1p                          12.17383
##                      max.Pplr.fctr.All.X..rcv.glmnet.Y
## PubDate.juliandate                           334.00000
## PubDate.last32.log1p                          12.21791
##                      min.Pplr.fctr.All.X..rcv.glmnet.N
## PubDate.juliandate                          244.000000
## PubDate.last32.log1p                          9.194211
##                      min.Pplr.fctr.All.X..rcv.glmnet.Y
## PubDate.juliandate                          244.000000
## PubDate.last32.log1p                          9.151227
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
## [1] "newobs Pplr.fctr.Final..rcv.glmnet Y: min < min of Train range: 10"
##      UniqueID Pplr.fctr.Final..rcv.glmnet PubDate.last32.log1p
## 6635     6635                           Y            10.171451
## 6636     6636                           Y            10.162461
## 7742     7742                           Y             9.087268
## 7743     7743                           Y             9.096836
## 7746     7746                           Y             9.060215
## 7748     7748                           Y             9.049115
## 7759     7759                           Y             9.069813
## 7938     7938                           Y             8.860641
## 7939     7939                           Y             8.835792
## 8217     8217                           Y            11.168659
##      PubDate.last8.log1p WordCount.log1p WordCount.root2
## 6635            7.156177        7.609862        44.91102
## 6636            6.891626        6.853299        30.75711
## 7742            7.708860        7.565275        43.92038
## 7743            7.748891        6.144186        21.56386
## 7746            7.610358        6.363028        24.06242
## 7748            7.536897        5.455321        15.26434
## 7759            7.826842        6.809039        30.08322
## 7938            7.494430        6.543912        26.34388
## 7939            7.101676        6.115892        21.26029
## 8217            9.528794        1.609438         2.00000
##                                        id      cor.y exclude.as.feat
## PubDate.last32.log1p PubDate.last32.log1p 0.02254251           FALSE
## PubDate.last8.log1p   PubDate.last8.log1p 0.05657458           FALSE
## WordCount.log1p           WordCount.log1p 0.25431963           FALSE
## WordCount.root2           WordCount.root2 0.29212068           FALSE
##                       cor.y.abs          cor.high.X freqRatio
## PubDate.last32.log1p 0.02254251                <NA>  1.000000
## PubDate.last8.log1p  0.05657458 PubDate.last4.log1p  1.166667
## WordCount.log1p      0.25431963     WordCount.root2  2.315789
## WordCount.root2      0.29212068                <NA>  2.315789
##                      percentUnique zeroVar   nzv is.cor.y.abs.low
## PubDate.last32.log1p      91.13595   FALSE FALSE            FALSE
## PubDate.last8.log1p       75.15309   FALSE FALSE            FALSE
## WordCount.log1p           24.15799   FALSE FALSE            FALSE
## WordCount.root2           24.15799   FALSE FALSE            FALSE
##                      interaction.feat shapiro.test.p.value rsp_var_raw
## PubDate.last32.log1p               NA         1.048418e-44       FALSE
## PubDate.last8.log1p                NA         1.021260e-42       FALSE
## WordCount.log1p                    NA         1.576866e-49       FALSE
## WordCount.root2                    NA         4.556481e-30       FALSE
##                      id_var rsp_var       max      min max.Pplr.fctr.N
## PubDate.last32.log1p     NA      NA  12.32341 8.835792       12.305461
## PubDate.last8.log1p      NA      NA  11.62246 6.666957       11.622461
## WordCount.log1p          NA      NA   9.29771 0.000000        8.819665
## WordCount.root2          NA      NA 104.46052 0.000000       82.249620
##                      max.Pplr.fctr.Y min.Pplr.fctr.N min.Pplr.fctr.Y
## PubDate.last32.log1p        12.21791        9.125654        9.176990
## PubDate.last8.log1p         11.44336        6.666957        7.187657
## WordCount.log1p              9.29771        0.000000        1.945910
## WordCount.root2            104.46052        0.000000        2.449490
##                      max.Pplr.fctr.All.X..rcv.glmnet.N
## PubDate.last32.log1p                         12.173834
## PubDate.last8.log1p                          11.401502
## WordCount.log1p                               7.059618
## WordCount.root2                              34.102786
##                      max.Pplr.fctr.All.X..rcv.glmnet.Y
## PubDate.last32.log1p                         12.217912
## PubDate.last8.log1p                          11.622461
## WordCount.log1p                               9.140883
## WordCount.root2                              96.581572
##                      min.Pplr.fctr.All.X..rcv.glmnet.N
## PubDate.last32.log1p                          9.194211
## PubDate.last8.log1p                           6.934397
## WordCount.log1p                               0.000000
## WordCount.root2                               0.000000
##                      min.Pplr.fctr.All.X..rcv.glmnet.Y
## PubDate.last32.log1p                          9.151227
## PubDate.last8.log1p                           7.019297
## WordCount.log1p                               0.000000
## WordCount.root2                               0.000000
##                      max.Pplr.fctr.Final..rcv.glmnet.N
## PubDate.last32.log1p                         12.323407
## PubDate.last8.log1p                          11.279555
## WordCount.log1p                               7.973844
## WordCount.root2                              53.879495
##                      max.Pplr.fctr.Final..rcv.glmnet.Y
## PubDate.last32.log1p                         12.305461
## PubDate.last8.log1p                          11.332279
## WordCount.log1p                               8.692322
## WordCount.root2                              77.175126
##                      min.Pplr.fctr.Final..rcv.glmnet.N
## PubDate.last32.log1p                          8.862908
## PubDate.last8.log1p                           7.061334
## WordCount.log1p                               0.000000
## WordCount.root2                               0.000000
##                      min.Pplr.fctr.Final..rcv.glmnet.Y
## PubDate.last32.log1p                          8.835792
## PubDate.last8.log1p                           6.891626
## WordCount.log1p                               1.609438
## WordCount.root2                               2.000000
## [1] "newobs Pplr.fctr.Final..rcv.glmnet Y: max > max of Train range: 681"
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
## 6549     6549                           Y               0.0103067569
## 6568     6568                           Y               0.0058403094
## 6859     6859                           Y              -0.0003328293
## 8074     8074                           Y              -0.0039277748
## 8137     8137                           Y              -0.0057434038
## 8179     8179                           Y              -0.0037098993
##      PubDate.day.minutes.poly.3 PubDate.day.minutes.poly.4
## 6549               -0.009159941              -0.0099955355
## 6568               -0.008760997              -0.0002028669
## 6859               -0.001377328               0.0092461622
## 8074                0.004086564               0.0087598734
## 8137                0.006605125               0.0066524422
## 8179                0.003766018               0.0089316619
##      PubDate.day.minutes.poly.5 PubDate.juliandate PubDate.last32.log1p
## 6549               0.0006942245                335             9.471242
## 6568               0.0095286285                335             9.661734
## 6859               0.0021075916                338             9.886443
## 8074              -0.0065295954                356            10.655847
## 8137              -0.0098074102                357            10.956440
## 8179              -0.0060564299                358            11.239304
##      WordCount.nexp
## 6549   0.000000e+00
## 6568  5.688906e-247
## 6859  7.149792e-142
## 8074   0.000000e+00
## 8137   0.000000e+00
## 8179   0.000000e+00
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
##                                                    id       cor.y
## PubDate.day.minutes.poly.1 PubDate.day.minutes.poly.1  0.15675348
## PubDate.day.minutes.poly.3 PubDate.day.minutes.poly.3  0.02798355
## PubDate.day.minutes.poly.4 PubDate.day.minutes.poly.4  0.07394139
## PubDate.day.minutes.poly.5 PubDate.day.minutes.poly.5 -0.05592923
## PubDate.juliandate                 PubDate.juliandate  0.01436107
## PubDate.last32.log1p             PubDate.last32.log1p  0.02254251
## WordCount.nexp                         WordCount.nexp -0.05320840
##                            exclude.as.feat  cor.y.abs         cor.high.X
## PubDate.day.minutes.poly.1           FALSE 0.15675348               <NA>
## PubDate.day.minutes.poly.3           FALSE 0.02798355               <NA>
## PubDate.day.minutes.poly.4           FALSE 0.07394139               <NA>
## PubDate.day.minutes.poly.5           FALSE 0.05592923               <NA>
## PubDate.juliandate                   FALSE 0.01436107 PubDate.month.fctr
## PubDate.last32.log1p                 FALSE 0.02254251               <NA>
## WordCount.nexp                       FALSE 0.05320840               <NA>
##                            freqRatio percentUnique zeroVar   nzv
## PubDate.day.minutes.poly.1   1.22549     18.080220   FALSE FALSE
## PubDate.day.minutes.poly.3   1.22549     18.080220   FALSE FALSE
## PubDate.day.minutes.poly.4   1.22549     18.080220   FALSE FALSE
## PubDate.day.minutes.poly.5   1.22549     18.080220   FALSE FALSE
## PubDate.juliandate           1.03252      1.393141   FALSE FALSE
## PubDate.last32.log1p         1.00000     91.135946   FALSE FALSE
## WordCount.nexp              17.76136     11.328843   FALSE FALSE
##                            is.cor.y.abs.low interaction.feat
## PubDate.day.minutes.poly.1            FALSE               NA
## PubDate.day.minutes.poly.3            FALSE               NA
## PubDate.day.minutes.poly.4            FALSE               NA
## PubDate.day.minutes.poly.5            FALSE               NA
## PubDate.juliandate                    FALSE               NA
## PubDate.last32.log1p                  FALSE               NA
## WordCount.nexp                        FALSE               NA
##                            shapiro.test.p.value rsp_var_raw id_var rsp_var
## PubDate.day.minutes.poly.1         1.590362e-18       FALSE     NA      NA
## PubDate.day.minutes.poly.3         9.822405e-52       FALSE     NA      NA
## PubDate.day.minutes.poly.4         1.523136e-47       FALSE     NA      NA
## PubDate.day.minutes.poly.5         1.157500e-41       FALSE     NA      NA
## PubDate.juliandate                 1.389406e-35       FALSE     NA      NA
## PubDate.last32.log1p               1.048418e-44       FALSE     NA      NA
## WordCount.nexp                     9.108805e-94       FALSE     NA      NA
##                                     max          min max.Pplr.fctr.N
## PubDate.day.minutes.poly.1   0.02475916  -0.02749464      0.02472285
## PubDate.day.minutes.poly.3   0.05215301  -0.04512497      0.05182988
## PubDate.day.minutes.poly.4   0.06677441  -0.01832740      0.06610094
## PubDate.day.minutes.poly.5   0.08471756  -0.02450918      0.08344228
## PubDate.juliandate         365.00000000 244.00000000    334.00000000
## PubDate.last32.log1p        12.32340669   8.83579237     12.30546086
## WordCount.nexp               1.00000000   0.00000000      1.00000000
##                            max.Pplr.fctr.Y min.Pplr.fctr.N min.Pplr.fctr.Y
## PubDate.day.minutes.poly.1    2.446866e-02     -0.02749464     -0.02745833
## PubDate.day.minutes.poly.3    4.959703e-02     -0.04512497     -0.04482024
## PubDate.day.minutes.poly.4    6.149053e-02     -0.01832740     -0.01821959
## PubDate.day.minutes.poly.5    7.481472e-02     -0.02450918     -0.02362780
## PubDate.juliandate            3.340000e+02    244.00000000    244.00000000
## PubDate.last32.log1p          1.221791e+01      9.12565356      9.17699039
## WordCount.nexp                2.478752e-03      0.00000000      0.00000000
##                            max.Pplr.fctr.All.X..rcv.glmnet.N
## PubDate.day.minutes.poly.1                        0.02450498
## PubDate.day.minutes.poly.3                        0.04991290
## PubDate.day.minutes.poly.4                        0.06213811
## PubDate.day.minutes.poly.5                        0.07601554
## PubDate.juliandate                              332.00000000
## PubDate.last32.log1p                             12.17383350
## WordCount.nexp                                    1.00000000
##                            max.Pplr.fctr.All.X..rcv.glmnet.Y
## PubDate.day.minutes.poly.1                        0.02472285
## PubDate.day.minutes.poly.3                        0.05182988
## PubDate.day.minutes.poly.4                        0.06610094
## PubDate.day.minutes.poly.5                        0.08344228
## PubDate.juliandate                              334.00000000
## PubDate.last32.log1p                             12.21791228
## WordCount.nexp                                    1.00000000
##                            min.Pplr.fctr.All.X..rcv.glmnet.N
## PubDate.day.minutes.poly.1                       -0.02749464
## PubDate.day.minutes.poly.3                       -0.04512497
## PubDate.day.minutes.poly.4                       -0.01832685
## PubDate.day.minutes.poly.5                       -0.02450918
## PubDate.juliandate                              244.00000000
## PubDate.last32.log1p                              9.19421099
## WordCount.nexp                                    0.00000000
##                            min.Pplr.fctr.All.X..rcv.glmnet.Y
## PubDate.day.minutes.poly.1                       -0.02745833
## PubDate.day.minutes.poly.3                       -0.04482024
## PubDate.day.minutes.poly.4                       -0.01824013
## PubDate.day.minutes.poly.5                       -0.02362780
## PubDate.juliandate                              244.00000000
## PubDate.last32.log1p                              9.15122711
## WordCount.nexp                                    0.00000000
##                            max.Pplr.fctr.Final..rcv.glmnet.N
## PubDate.day.minutes.poly.1                        0.02468654
## PubDate.day.minutes.poly.3                        0.05150779
## PubDate.day.minutes.poly.4                        0.06543120
## PubDate.day.minutes.poly.5                        0.08217780
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
##            MFO###myMFO_classfr      Random###myrandom_classfr 
##                              0                              0 
##     Max.cor.Y.rcv.1X1###glmnet Max.cor.Y.rcv.1X1.cp.0###rpart 
##                              0                              0
```

```
##                                 max.Accuracy.OOB max.AUCROCR.OOB
## Max.cor.Y##rcv#rpart                   0.8200231       0.5892132
## Max.cor.Y.Time.Poly##rcv#glmnet        0.7754630       0.8049472
## Interact.High.cor.Y##rcv#glmnet        0.7743056       0.7992251
## Max.cor.Y.rcv.1X1.cp.0###rpart         0.7673611       0.7773858
## Max.cor.Y.rcv.1X1###glmnet             0.7604167       0.8116126
## Max.cor.Y.rcv.5X3##rcv#glmnet          0.7604167       0.8114863
## Max.cor.Y.rcv.5X1##rcv#glmnet          0.7604167       0.8114863
## Max.cor.Y.rcv.5X5##rcv#glmnet          0.7604167       0.8114863
## Max.cor.Y.rcv.3X1##rcv#glmnet          0.7575231       0.8067975
## Max.cor.Y.rcv.3X3##rcv#glmnet          0.7575231       0.8067975
## Max.cor.Y.rcv.3X5##rcv#glmnet          0.7575231       0.8067975
## Max.cor.Y.Time.Lag##rcv#glmnet         0.6464120       0.8118709
## Low.cor.X##rcv#glmnet                  0.6377315       0.8150026
## All.X##rcv#glmnet                      0.6261574       0.8157921
## MFO###myMFO_classfr                    0.1331019       0.5000000
## Random###myrandom_classfr              0.1331019       0.4857956
## Final##rcv#glmnet                             NA              NA
##                                 max.AUCpROC.OOB max.Accuracy.fit
## Max.cor.Y##rcv#rpart                  0.5870523        0.9296422
## Max.cor.Y.Time.Poly##rcv#glmnet       0.5965780        0.9322790
## Interact.High.cor.Y##rcv#glmnet       0.5957392        0.9315850
## Max.cor.Y.rcv.1X1.cp.0###rpart        0.6174697        0.9381765
## Max.cor.Y.rcv.1X1###glmnet            0.5962443        0.9329725
## Max.cor.Y.rcv.5X3##rcv#glmnet         0.5962443        0.9333905
## Max.cor.Y.rcv.5X1##rcv#glmnet         0.5962443        0.9331818
## Max.cor.Y.rcv.5X5##rcv#glmnet         0.5962443        0.9331816
## Max.cor.Y.rcv.3X1##rcv#glmnet         0.5962443        0.9335973
## Max.cor.Y.rcv.3X3##rcv#glmnet         0.5962443        0.9333193
## Max.cor.Y.rcv.3X5##rcv#glmnet         0.5962443        0.9332218
## Max.cor.Y.Time.Lag##rcv#glmnet        0.5928717        0.9281861
## Low.cor.X##rcv#glmnet                 0.5876850        0.9261046
## All.X##rcv#glmnet                     0.5876850        0.9263124
## MFO###myMFO_classfr                   0.5000000        0.1796420
## Random###myrandom_classfr             0.5125675        0.1796420
## Final##rcv#glmnet                            NA        0.9063073
##                                 opt.prob.threshold.fit
## Max.cor.Y##rcv#rpart                               0.6
## Max.cor.Y.Time.Poly##rcv#glmnet                    0.4
## Interact.High.cor.Y##rcv#glmnet                    0.4
## Max.cor.Y.rcv.1X1.cp.0###rpart                     0.4
## Max.cor.Y.rcv.1X1###glmnet                         0.5
## Max.cor.Y.rcv.5X3##rcv#glmnet                      0.5
## Max.cor.Y.rcv.5X1##rcv#glmnet                      0.5
## Max.cor.Y.rcv.5X5##rcv#glmnet                      0.5
## Max.cor.Y.rcv.3X1##rcv#glmnet                      0.4
## Max.cor.Y.rcv.3X3##rcv#glmnet                      0.4
## Max.cor.Y.rcv.3X5##rcv#glmnet                      0.4
## Max.cor.Y.Time.Lag##rcv#glmnet                     0.3
## Low.cor.X##rcv#glmnet                              0.3
## All.X##rcv#glmnet                                  0.3
## MFO###myMFO_classfr                                0.1
## Random###myrandom_classfr                          0.1
## Final##rcv#glmnet                                  0.2
##                                 opt.prob.threshold.OOB
## Max.cor.Y##rcv#rpart                               0.6
## Max.cor.Y.Time.Poly##rcv#glmnet                    0.1
## Interact.High.cor.Y##rcv#glmnet                    0.1
## Max.cor.Y.rcv.1X1.cp.0###rpart                     0.1
## Max.cor.Y.rcv.1X1###glmnet                         0.1
## Max.cor.Y.rcv.5X3##rcv#glmnet                      0.1
## Max.cor.Y.rcv.5X1##rcv#glmnet                      0.1
## Max.cor.Y.rcv.5X5##rcv#glmnet                      0.1
## Max.cor.Y.rcv.3X1##rcv#glmnet                      0.1
## Max.cor.Y.rcv.3X3##rcv#glmnet                      0.1
## Max.cor.Y.rcv.3X5##rcv#glmnet                      0.1
## Max.cor.Y.Time.Lag##rcv#glmnet                     0.1
## Low.cor.X##rcv#glmnet                              0.1
## All.X##rcv#glmnet                                  0.1
## MFO###myMFO_classfr                                0.1
## Random###myrandom_classfr                          0.1
## Final##rcv#glmnet                                   NA
```

```
## [1] "All.X##rcv#glmnet OOB confusion matrix & accuracy: "
##          Prediction
## Reference   N   Y
##         N 874 624
##         Y  22 208
##                         .freqRatio.Fit .freqRatio.OOB .freqRatio.Tst
## OpEd#Opnn#                 0.090965862   0.0515046296    0.087700535
## #Opnn#ThPblcEdtr           0.003330558   0.0023148148    0.005347594
## Styls#U.S.#                0.026436303   0.0289351852    0.032620321
## Bsnss#Crsswrds/Gms#        0.021856786   0.0104166667    0.022459893
## Scnc#Hlth#                 0.030807660   0.0277777778    0.030481283
## Bsnss#Tchnlgy#             0.044338052   0.0729166667    0.060962567
## ##                         0.190049958   0.2146990741    0.182887701
## Bsnss#BsnssDy#Dlbk         0.130932556   0.1869212963    0.162566845
## Mtr#N.Y./Rgn#              0.026644463   0.0405092593    0.035828877
## Cltr#Arts#                 0.101998335   0.1070601852    0.093048128
## #Opnn#RmFrDbt              0.008742714   0.0115740741    0.010695187
## Styls##Fshn                0.021648626   0.0086805556    0.008021390
## Bsnss#BsnssDy#SmllBsnss    0.020815987   0.0231481481    0.021925134
## myOthr                     0.006869276   0.0028935185    0.002673797
## Trvl#Trvl#                 0.017277269   0.0196759259    0.018716578
## Cltr##                              NA   0.0005787037    0.037433155
## Frgn#Wrld#AsPcfc           0.031223980   0.0306712963    0.029946524
## #Mltmd#                    0.019150708   0.0283564815    0.027807487
## TStyl##                    0.129683597   0.0584490741    0.056149733
## #U.S.#Edctn                0.050582848   0.0474537037    0.047593583
## Frgn#Wrld#                 0.026644463   0.0254629630    0.025133690
##                         .n.Fit .n.New.N .n.New.Y .n.OOB .n.Trn.N .n.Trn.Y
## OpEd#Opnn#                 437       NA      164     89      117      409
## #Opnn#ThPblcEdtr            16       NA       10      4        4       16
## Styls#U.S.#                127       NA       61     50       77      100
## Bsnss#Crsswrds/Gms#        105       NA       42     18       20      103
## Scnc#Hlth#                 148       NA       57     48       74      122
## Bsnss#Tchnlgy#             213       34       80    126      288       51
## ##                         913      262       80    371     1169      115
## Bsnss#BsnssDy#Dlbk         629      201      103    323      864       88
## Mtr#N.Y./Rgn#              128       36       31     70      181       17
## Cltr#Arts#                 490      157       17    185      625       50
## #Opnn#RmFrDbt               42       20       NA     20       61        1
## Styls##Fshn                104       15       NA     15      118        1
## Bsnss#BsnssDy#SmllBsnss    100       36        5     40      135        5
## myOthr                      33        5       NA      5       38       NA
## Trvl#Trvl#                  83       35       NA     34      116        1
## Cltr##                      NA       47       23      1        1       NA
## Frgn#Wrld#AsPcfc           150       51        5     53      200        3
## #Mltmd#                     92       52       NA     49      139        2
## TStyl##                    623      104        1    101      715        9
## #U.S.#Edctn                243       87        2     82      325       NA
## Frgn#Wrld#                 128       47       NA     44      172       NA
##                         .n.Tst .n.fit .n.new .n.trn err.abs.OOB.mean
## OpEd#Opnn#                 164    437    164    526       0.52330093
## #Opnn#ThPblcEdtr            10     16     10     20       0.50201689
## Styls#U.S.#                 61    127     61    177       0.47265683
## Bsnss#Crsswrds/Gms#         42    105     42    123       0.46617390
## Scnc#Hlth#                  57    148     57    196       0.46081076
## Bsnss#Tchnlgy#             114    213    114    339       0.23596749
## ##                         342    913    342   1284       0.21028224
## Bsnss#BsnssDy#Dlbk         304    629    304    952       0.20486485
## Mtr#N.Y./Rgn#               67    128     67    198       0.19249438
## Cltr#Arts#                 174    490    174    675       0.18753758
## #Opnn#RmFrDbt               20     42     20     62       0.18712633
## Styls##Fshn                 15    104     15    119       0.14376823
## Bsnss#BsnssDy#SmllBsnss     41    100     41    140       0.14054271
## myOthr                       5     33      5     38       0.11351975
## Trvl#Trvl#                  35     83     35    117       0.10662343
## Cltr##                      70     NA     70      1       0.10322006
## Frgn#Wrld#AsPcfc            56    150     56    203       0.10088812
## #Mltmd#                     52     92     52    141       0.09869284
## TStyl##                    105    623    105    724       0.09409437
## #U.S.#Edctn                 89    243     89    325       0.07342928
## Frgn#Wrld#                  47    128     47    172       0.07272040
##                         err.abs.fit.mean err.abs.new.mean err.abs.trn.mean
## OpEd#Opnn#                    0.38883472               NA       0.33295283
## #Opnn#ThPblcEdtr              0.45394245               NA       0.32096412
## Styls#U.S.#                   0.49068821               NA       0.45281920
## Bsnss#Crsswrds/Gms#           0.35845239               NA       0.23739811
## Scnc#Hlth#                    0.45418853               NA       0.38566104
## Bsnss#Tchnlgy#                0.22005479               NA       0.23159760
## ##                            0.14540647               NA       0.13629887
## Bsnss#BsnssDy#Dlbk            0.15383991               NA       0.15219558
## Mtr#N.Y./Rgn#                 0.15574134               NA       0.15274373
## Cltr#Arts#                    0.12292240               NA       0.10912086
## #Opnn#RmFrDbt                 0.15566721               NA       0.05718416
## Styls##Fshn                   0.08609432               NA       0.03338503
## Bsnss#BsnssDy#SmllBsnss       0.13032416               NA       0.08239213
## myOthr                        0.11282519               NA       0.02825053
## Trvl#Trvl#                    0.08193916               NA       0.03226070
## Cltr##                                NA               NA       0.07220822
## Frgn#Wrld#AsPcfc              0.10382572               NA       0.04236734
## #Mltmd#                       0.09099175               NA       0.04168488
## TStyl##                       0.06990527               NA       0.02853194
## #U.S.#Edctn                   0.06377772               NA       0.01167100
## Frgn#Wrld#                    0.06994962               NA       0.01407013
##                         err.abs.OOB.sum err.abs.fit.sum err.abs.new.sum
## OpEd#Opnn#                   46.5737830      169.920774              NA
## #Opnn#ThPblcEdtr              2.0080675        7.263079              NA
## Styls#U.S.#                  23.6328413       62.317403              NA
## Bsnss#Crsswrds/Gms#           8.3911301       37.637501              NA
## Scnc#Hlth#                   22.1189167       67.219902              NA
## Bsnss#Tchnlgy#               29.7319043       46.871671              NA
## ##                           78.0147094      132.756111              NA
## Bsnss#BsnssDy#Dlbk           66.1713471       96.765300              NA
## Mtr#N.Y./Rgn#                13.4746063       19.934891              NA
## Cltr#Arts#                   34.6944525       60.231977              NA
## #Opnn#RmFrDbt                 3.7425266        6.538023              NA
## Styls##Fshn                   2.1565235        8.953809              NA
## Bsnss#BsnssDy#SmllBsnss       5.6217083       13.032416              NA
## myOthr                        0.5675988        3.723231              NA
## Trvl#Trvl#                    3.6251966        6.800950              NA
## Cltr##                        0.1032201              NA              NA
## Frgn#Wrld#AsPcfc              5.3470704       15.573858              NA
## #Mltmd#                       4.8359493        8.371241              NA
## TStyl##                       9.5035315       43.550981              NA
## #U.S.#Edctn                   6.0212009       15.497987              NA
## Frgn#Wrld#                    3.1996976        8.953551              NA
##                         err.abs.trn.sum
## OpEd#Opnn#                 175.13319010
## #Opnn#ThPblcEdtr             6.41928249
## Styls#U.S.#                 80.14899777
## Bsnss#Crsswrds/Gms#         29.19996741
## Scnc#Hlth#                  75.58956367
## Bsnss#Tchnlgy#              78.51158544
## ##                         175.00774839
## Bsnss#BsnssDy#Dlbk         144.89018795
## Mtr#N.Y./Rgn#               30.24325893
## Cltr#Arts#                  73.65658273
## #Opnn#RmFrDbt                3.54541784
## Styls##Fshn                  3.97281808
## Bsnss#BsnssDy#SmllBsnss     11.53489794
## myOthr                       1.07352023
## Trvl#Trvl#                   3.77450215
## Cltr##                       0.07220822
## Frgn#Wrld#AsPcfc             8.60057057
## #Mltmd#                      5.87756799
## TStyl##                     20.65712219
## #U.S.#Edctn                  3.79307552
## Frgn#Wrld#                   2.42006229
##   .freqRatio.Fit   .freqRatio.OOB   .freqRatio.Tst           .n.Fit 
##               NA         1.000000         1.000000               NA 
##         .n.New.N         .n.New.Y           .n.OOB         .n.Trn.N 
##               NA               NA      1728.000000      5439.000000 
##         .n.Trn.Y           .n.Tst           .n.fit           .n.new 
##               NA      1870.000000               NA      1870.000000 
##           .n.trn err.abs.OOB.mean err.abs.fit.mean err.abs.new.mean 
##      6532.000000         4.690731               NA               NA 
## err.abs.trn.mean  err.abs.OOB.sum  err.abs.fit.sum  err.abs.new.sum 
##         2.955758       369.535982               NA               NA 
##  err.abs.trn.sum 
##       934.122128
```

```
##                                         All.X__rcv_glmnet.imp
## PubDate.day.minutes.poly.1                         100.000000
## PubDate.day.minutes.poly.4                          46.118345
## NDSSName.my.fctrOpEd#Opnn#                          29.608813
## NDSSName.my.fctrBsnss#Crsswrds/Gms#                 27.135657
## NDSSName.my.fctrScnc#Hlth#                          24.964890
## NDSSName.my.fctr#Opnn#ThPblcEdtr                    24.849247
## NDSSName.my.fctrStyls#U.S.#                         23.449567
## PubDate.day.minutes.poly.2                          22.879214
## WordCount.log1p                                      7.599061
## PubDate.wkend                                        7.555716
## PubDate.hour.fctr(15.3,23]                           6.568901
## WordCount.root2                                      6.405483
## PubDate.last4.log1p                                  6.403279
## PubDate.last2.log1p                                  6.307833
## PubDate.last8.log1p                                  6.231915
## PubDate.day.minutes.poly.3                           6.181073
## NDSSName.my.fctrBsnss#Tchnlgy#                       6.181073
## PubDate.hour.fctr(7.67,15.3]                         6.181073
## PubDate.wkday.fctr1                                  6.181073
## PubDate.last16.log1p                                 6.181073
## PubDate.date.fctr(7,13]                              6.181073
## PubDate.date.fctr(25,31]                             6.181073
## .rnorm                                               6.181073
## NDSSName.my.fctrCltr##                               6.181073
## NDSSName.my.fctrMtr#N.Y./Rgn#                        6.181073
## PubDate.date.fctr(19,25]                             6.181073
## PubDate.day.minutes.poly.5                           6.181073
## PubDate.last32.log1p                                 6.181073
## PubDate.minute.fctr(14.8,29.5]                       6.181073
## PubDate.minute.fctr(44.2,59.1]                       6.181073
## PubDate.month.fctr10                                 6.181073
## PubDate.month.fctr12                                 6.181073
## PubDate.second.fctr(14.8,29.5]                       6.181073
## PubDate.second.fctr(29.5,44.2]                       6.181073
## PubDate.wkday.fctr3                                  6.181073
## PubDate.wkday.fctr4                                  6.181073
## WordCount.nexp                                       6.181073
## PubDate.juliandate                                   6.181073
## PubDate.wkday.fctr2                                  6.181073
## PubDate.date.fctr(13,19]                             6.181073
## PubDate.month.fctr11                                 6.181073
## PubDate.second.fctr(44.2,59.1]                       6.181073
## PubDate.wkday.fctr5                                  6.181073
## PubDate.minute.fctr(29.5,44.2]                       6.181073
## PubDate.wkday.fctr6                                  6.181073
## NDSSName.my.fctrTrvl#Trvl#                           4.901316
## NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss              4.791503
## NDSSName.my.fctrBsnss#BsnssDy#Dlbk                   4.575817
## NDSSName.my.fctrCltr#Arts#                           4.494666
##                                         Final__rcv_glmnet.imp
## PubDate.day.minutes.poly.1                           78.09046
## PubDate.day.minutes.poly.4                           21.36267
## NDSSName.my.fctrOpEd#Opnn#                           29.73570
## NDSSName.my.fctrBsnss#Crsswrds/Gms#                  27.53840
## NDSSName.my.fctrScnc#Hlth#                           25.24514
## NDSSName.my.fctr#Opnn#ThPblcEdtr                     27.29376
## NDSSName.my.fctrStyls#U.S.#                          24.24782
## PubDate.day.minutes.poly.2                          100.00000
## WordCount.log1p                                      15.19157
## PubDate.wkend                                        15.37916
## PubDate.hour.fctr(15.3,23]                           13.49353
## WordCount.root2                                      13.72215
## PubDate.last4.log1p                                  13.49353
## PubDate.last2.log1p                                  13.55706
## PubDate.last8.log1p                                  13.49353
## PubDate.day.minutes.poly.3                           28.57451
## NDSSName.my.fctrBsnss#Tchnlgy#                       15.43897
## PubDate.hour.fctr(7.67,15.3]                         14.61230
## PubDate.wkday.fctr1                                  13.91542
## PubDate.last16.log1p                                 13.80029
## PubDate.date.fctr(7,13]                              13.67140
## PubDate.date.fctr(25,31]                             13.63428
## .rnorm                                               13.49353
## NDSSName.my.fctrCltr##                               13.49353
## NDSSName.my.fctrMtr#N.Y./Rgn#                        13.49353
## PubDate.date.fctr(19,25]                             13.49353
## PubDate.day.minutes.poly.5                           13.49353
## PubDate.last32.log1p                                 13.49353
## PubDate.minute.fctr(14.8,29.5]                       13.49353
## PubDate.minute.fctr(44.2,59.1]                       13.49353
## PubDate.month.fctr10                                 13.49353
## PubDate.month.fctr12                                 13.49353
## PubDate.second.fctr(14.8,29.5]                       13.49353
## PubDate.second.fctr(29.5,44.2]                       13.49353
## PubDate.wkday.fctr3                                  13.49353
## PubDate.wkday.fctr4                                  13.49353
## WordCount.nexp                                       13.49353
## PubDate.juliandate                                   13.49350
## PubDate.wkday.fctr2                                  13.31814
## PubDate.date.fctr(13,19]                             13.26769
## PubDate.month.fctr11                                 13.16393
## PubDate.second.fctr(44.2,59.1]                       13.11401
## PubDate.wkday.fctr5                                  12.82054
## PubDate.minute.fctr(29.5,44.2]                       12.78565
## PubDate.wkday.fctr6                                  12.51262
## NDSSName.my.fctrTrvl#Trvl#                           10.36195
## NDSSName.my.fctrBsnss#BsnssDy#SmllBsnss              10.18924
## NDSSName.my.fctrBsnss#BsnssDy#Dlbk                   12.35507
## NDSSName.my.fctrCltr#Arts#                           12.63313
```

```
## [1] "glbObsNew prediction stats:"
```

```
## 
##    N    Y 
## 1189  681
```

```
##                   label step_major step_minor label_minor     bgn     end
## 16     predict.data.new          8          0           0 358.247 379.494
## 17 display.session.info          9          0           0 379.495      NA
##    elapsed
## 16  21.248
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
## 10              fit.models          6          0           0  83.208
## 14       fit.data.training          7          0           0 282.564
## 5         extract.features          3          0           0  24.453
## 11              fit.models          6          1           1 231.997
## 16        predict.data.new          8          0           0 358.247
## 9          select.features          5          0           0  62.061
## 12              fit.models          6          2           2 261.729
## 1              import.data          1          0           0   9.172
## 15       fit.data.training          7          1           1 347.894
## 13              fit.models          6          3           3 276.543
## 2             inspect.data          2          0           0  20.292
## 8  partition.data.training          4          0           0  60.746
## 6      manage.missing.data          3          1           1  59.576
## 3               scrub.data          2          1           1  23.327
## 4           transform.data          2          2           2  24.356
## 7             cluster.data          3          2           2  60.677
##        end elapsed duration
## 10 231.997 148.789  148.789
## 14 347.894  65.330   65.330
## 5   59.576  35.123   35.123
## 11 261.728  29.731   29.731
## 16 379.494  21.248   21.247
## 9   83.208  21.147   21.147
## 12 276.542  14.813   14.813
## 1   20.291  11.119   11.119
## 15 358.247  10.353   10.353
## 13 282.563   6.020    6.020
## 2   23.326   3.034    3.034
## 8   62.060   1.314    1.314
## 6   60.676   1.101    1.100
## 3   24.356   1.029    1.029
## 4   24.452   0.096    0.096
## 7   60.745   0.068    0.068
## [1] "Total Elapsed Time: 379.494 secs"
```

![](NYTBlogs3_feat_PubDate_files/figure-html/display.session.info-1.png) 

```
##                                   label step_major step_minor
## 4        fit.models_0_Max.cor.Y.rcv.*X*          1          3
## 5 fit.models_0_Max.cor.Y[rcv.1X1.cp.0|]          1          4
## 9                fit.models_0_Low.cor.X          1          8
## 8      fit.models_0_Interact.High.cor.Y          1          7
## 6      fit.models_0_Max.cor.Y.Time.Poly          1          5
## 7       fit.models_0_Max.cor.Y.Time.Lag          1          6
## 3                   fit.models_0_Random          1          2
## 2                      fit.models_0_MFO          1          1
## 1                      fit.models_0_bgn          1          0
##        label_minor     bgn     end elapsed duration
## 4           glmnet  91.766 169.327  77.561   77.561
## 5            rpart 169.327 183.585  14.258   14.258
## 9           glmnet 218.981 231.982  13.001   13.001
## 8           glmnet 206.243 218.980  12.738   12.737
## 6           glmnet 183.586 195.343  11.757   11.757
## 7           glmnet 195.344 206.242  10.899   10.898
## 3 myrandom_classfr  87.427  91.766   4.339    4.339
## 2    myMFO_classfr  84.356  87.427   3.071    3.071
## 1            setup  84.325  84.355   0.031    0.030
## [1] "Total Elapsed Time: 231.982 secs"
```

![](NYTBlogs3_feat_PubDate_files/figure-html/display.session.info-2.png) 
