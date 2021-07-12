##################################################################################################
###                                                                                              #
###        This Script is used to do age prediction based on DNA methylation data(450K)          #
###        Coefficients of the predictor are based on 13,566 training samples                    #
###        Two age predictors based on different methods(Elastic Net and BLUP) can be used       #
###        Qian Zhang 27-03-2018, Email: q.zhang@uq.edu.au                                       #
##################################################################################################

library(ggplot2)

############# for each probe, change to missing value to the mean value across all individuals #############
addna <- function(methy) {
        methy[is.na(methy)] <- mean(methy, na.rm = TRUE)
        return(methy)
}

age.predictor <- function(data) {
        ############# 3. get the coefficients of each probe from Elastic Net/BLUP method, !!!!WE HAVE TWO PREDICTORS!!!#############
        
        print("1. Loading predictors")
        encoef <- read.table(
                "code/DNAm-based-age-predictor/en.coef",
                stringsAsFactor = F,
                header = T
        ) 
        blupcoef <- read.table(
                "code/DNAm-based-age-predictor/blup.coef",
                stringsAsFactor = F,
                header = T
        ) 
        
        en_int <- encoef[1, 2]
        blup_int <- blupcoef[1, 2]
        
        encoef <- encoef[-1, ]
        blupcoef <- blupcoef[-1, ]
        
        rownames(encoef) <- encoef$probe
        rownames(blupcoef) <- blupcoef$probe
        
        
        ############# 2. data loading and QC ##################
        print("2. Data loading and QC")
        
        if (nrow(data) > ncol(data)) {
                print("I guess you are using Probe in the row, data will be transformed!!!")
                data <- t(data)
        }
        
        print("2.1 Replacing missing values with mean value")
        dataNona <- apply(data, 2, function(x)  addna(x))   #  replace the NA with mean value for each probe
        
        print("2.2 Removing probe that has NA across all individuals")
        dataNona <-
                dataNona[, colSums(is.na(dataNona)) != nrow(dataNona)]
        print(
                paste0(
                        ncol(data) - ncol(dataNona),
                        " probe(s) is(are) removed since it has (they have) NA across all individuals"
                )
        )
        
        ############### standardize the DNA methylation within each individual,
        # remove the mean and divided by the SD of each individual Probe * IND
        print("2.3 Standardizing")
        dataNona.norm <- apply(dataNona, 1, scale)
        rownames(dataNona.norm) <- colnames(dataNona)
        
        ############# 4. get common probes between predictors and data ##############
        print("3. Checking misssing probes")
        
        encomm <- intersect(rownames(encoef), rownames(dataNona.norm))
        blupcomm <- intersect(rownames(blupcoef), rownames(dataNona.norm))
        
        endiff <- nrow(encoef) - length(encomm)
        blupdiff <- nrow(blupcoef) - length(blupcomm)
        
        print(
                paste0(
                        endiff,
                        " probe(s) in Elastic Net predictor is(are) not in the data"
                )
        )
        print(
                paste0(
                        blupdiff,
                        " probe(s) in BLUP predictor is(are) not in the data"
                )
        )
        print("BLUP can perform better if the number of missing probes is too large!")
        
        ############# 5. extract the common probes and do age prediction ###############
        print("4. Predicting")
        
        encoef <- encoef[encomm, ]
        blupcoef <- blupcoef[blupcomm, ]
        encoef$coef %*% dataNona.norm[encomm, ] + en_int -> enpred
        blupcoef$coef %*% dataNona.norm[blupcomm, ] + blup_int -> blupred
        
        ############# 6. Save the predicted result ###########
        age.pred <- data.frame(
                "enpred" =  as.double(enpred),
                "blupred" = as.double(blupred)
        )
        ################ DONE ########################################
        
        print("Completed!!!")
        
        ############ This is the script for plotting figure and calculating prediction accuracy ##########
        library(ggplot2)
        library(reshape2)
        colnames(age.pred) <- c("Elastic_Net", "BLUP")
        # age <- melt(age, measure.vars = c("Elastic Net", "BLUP"))
        plot.pred <- ggplot(data = age.pred, aes(x = Elastic_Net, y = BLUP)) +
                geom_abline(intercept = 0, slope = 1) +
                geom_point() +
                xlab("Predicted Age - Elastic_Net") +
                ylab("Predicted Age - BLUP") +
                theme(
                        axis.title = element_text(size = 12, face = "bold"),
                        legend.text = element_text(size = 10, face = "bold"),
                        legend.title = element_text(face = "bold"),
                        axis.text = element_text(size = 10, face = "bold"),
                        strip.text.x = element_text(face = "bold")
                )
        rownames(age.pred) <- rownames(data)
        return(
                list("age.pred" = age.pred, "plot.pred" = plot.pred)
        )
}
