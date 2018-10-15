#' Extract clinical information from tcga format
#' @export
extract_tcga_clinical <- function(clin, ESAC = FALSE, STAD = FALSE) {
    ## age and gender patient.age_at_initial_pathologic_diagnosis
    if (ESAC) {
        age <- clin[, "patient.primary_pathology.age_at_initial_pathologic_diagnosis"]
    } else {
        age <- clin[, "patient.age_at_initial_pathologic_diagnosis"]
    }
    gender <- clin[, "patient.gender"]
    
    # chemo_resp
    chemo_resp <- clin[, "patient.follow_ups.follow_up.primary_therapy_outcome_success"]
    names(chemo_resp) <- toupper(clin[, "patient.bcr_patient_barcode"])
    chemo_resp[which(chemo_resp == "complete remission/response")] <- sub("complete remission/response", 
        "sen", chemo_resp[which(chemo_resp == "complete remission/response")])
    chemo_resp[which(chemo_resp == "partial remission/response")] <- sub("partial remission/response", 
        "p_sen", chemo_resp[which(chemo_resp == "partial remission/response")])
    chemo_resp[which(chemo_resp == "stable disease")] <- sub("stable disease", 
        "resist", chemo_resp[which(chemo_resp == "stable disease")])
    chemo_resp[which(chemo_resp == "progressive disease")] <- sub("progressive disease", 
        "resist", chemo_resp[which(chemo_resp == "progressive disease")])
    
    # chemo
    chemo <- clin[, "patient.drugs.drug.therapy_types.therapy_type"]
    names(chemo) <- toupper(clin[, "patient.bcr_patient_barcode"])
    chemo[is.na(chemo)] <- 0
    chemo <- ifelse(chemo == "chemotherapy", 1, 0)
    
    ## t-stage
    Tstage <- clin[, "patient.stage_event.tnm_categories.pathologic_categories.pathologic_t"]
    Tstage[grep("t4", Tstage)] <- 4
    Tstage[grep("t3", Tstage)] <- 3
    Tstage[grep("t2", Tstage)] <- 2
    Tstage[grep("t1", Tstage)] <- 1
    Tstage <- as.numeric(Tstage)
    
    ## lympho node
    Nstage <- clin[, "patient.stage_event.tnm_categories.pathologic_categories.pathologic_n"]
    LN <- !Nstage == "n0"
    LN[Nstage == "nx"] <- NA
    
    ## stage
    stage <- clin[, "patient.stage_event.pathologic_stage"]
    names(stage) <- toupper(clin[, "patient.bcr_patient_barcode"])
    pt_stage23 <- names(stage)[grep("ii", stage)]
    stage[grep("iv", stage)] <- 4
    stage[grep("iii", stage)] <- 3
    stage[grep("ii", stage)] <- 2
    stage[grep("i", stage)] <- 1
    s <- as.numeric(stage)
    names(s) <- names(stage)
    stage <- s
    
    ## relapse
    if (ESAC) {
        labels <- clin[, "patient.new_tumor_events.new_tumor_event_after_initial_treatment"]
    } else if (STAD) {
        labels <- clin[, "patient.follow_ups.follow_up.new_tumor_event_after_initial_treatment"]
    } else {
        labels <- clin[, "patient.follow_ups.follow_up.new_tumor_events.new_tumor_event_after_initial_treatment"]
    }
    names(labels) <- toupper(clin[, "patient.bcr_patient_barcode"])
    relapse <- (labels == "yes")
    names(relapse) <- names(labels)
    
    ## relapse time
    if (ESAC) {
        t <- c("patient.new_tumor_events.new_tumor_event.days_to_new_tumor_event_after_initial_treatment", 
            "patient.days_to_death", "patient.follow_ups.follow_up.days_to_last_followup")
    } else if (STAD) {
        t <- c("patient.follow_ups.follow_up.days_to_new_tumor_event_after_initial_treatment", 
            "patient.follow_ups.follow_up.days_to_death", "patient.follow_ups.follow_up.days_to_last_followup")
    } else {
        t <- c("patient.follow_ups.follow_up.new_tumor_events.new_tumor_event.days_to_new_tumor_event_after_initial_treatment", 
            "patient.follow_ups.follow_up.days_to_death", "patient.follow_ups.follow_up.days_to_last_followup")
    }
    tmp <- clin[, t]
    colnames(tmp) <- NULL
    time <- tmp[, 1]
    time[is.na(time)] <- as.vector(t(tmp[is.na(time), 2]))
    time[is.na(time)] <- as.vector(t(tmp[is.na(time), 3]))
    relapse.time <- as.numeric(time)
    names(relapse.time) <- toupper(clin[, "patient.bcr_patient_barcode"])
    
    ## os
    labels <- clin[, "patient.vital_status"]
    names(labels) <- toupper(clin[, "patient.bcr_patient_barcode"])
    death <- (labels == "dead")
    names(death) <- names(labels)
    
    ## os time
    if (ESAC) {
        t <- c("patient.days_to_death", "patient.follow_ups.follow_up.days_to_last_followup")
    } else {
        t <- c("patient.follow_ups.follow_up.days_to_death", "patient.follow_ups.follow_up.days_to_last_followup")
    }
    tmp <- clin[, t]
    colnames(tmp) <- NULL
    time <- tmp[, 1]
    time[is.na(time)] <- as.vector(t(tmp[is.na(time), 2]))
    death.time <- as.numeric(time)
    names(death.time) <- toupper(clin[, "patient.bcr_patient_barcode"])
    
    # DFS
    DFS <- relapse | death
    DFS[is.na(DFS) & death] <- T
    
    if (ESAC) {
        type <- factor(clin[, "patient.biospecimen_cqcf.histological_type"] == 
            "esophagus adenocarcinoma  nos", levels = c(T, F), labels = c("EAC", 
            "ESCC"))
        location <- clin[, "patient.primary_pathology.esophageal_tumor_cental_location"]
        data.frame(age, gender, type, location, stage, LN, Tstage, DFS, 
            relapse, relapse.time, death, death.time)
    } else {
        data.frame(age, gender, stage, LN, Tstage, DFS, relapse, relapse.time, 
            death, death.time, chemo, chemo_resp)
    }
}
