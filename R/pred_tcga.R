#' Extract clinical information from tcga format
#' @export
extract_tcga_clinical <- function(clin) {
  ## Stage

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

  ## Relapse

  labels <- clin[, "patient.follow_ups.follow_up.new_tumor_events.new_tumor_event_after_initial_treatment"]
  names(labels) <- toupper(clin[, "patient.bcr_patient_barcode"])

  relapse <- (labels == "yes")
  names(relapse) <- names(labels)

  ## Time
  t <- c("patient.follow_ups.follow_up.new_tumor_events.new_tumor_event.days_to_new_tumor_event_after_initial_treatment",
         "patient.follow_ups.follow_up.days_to_death",
         "patient.follow_ups.follow_up.days_to_last_followup"
  )
  tmp <- clin[, t]
  colnames(tmp) <- NULL

  time <- tmp[, 1]
  time[is.na(time)] <- as.vector(t(tmp[is.na(time), 2]))
  time[is.na(time)] <- as.vector(t(tmp[is.na(time), 3]))
  time <- as.numeric(time)
  names(time) <- toupper(clin[, "patient.bcr_patient_barcode"])

  data.frame(stage, relapse, time)
}
