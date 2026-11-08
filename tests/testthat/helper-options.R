# Most fixtures fit the longitudinal sub-model on the subjects with an event,
# which makes predictRisk() warn once per session (check_complete_case_fit());
# test-fitImputedBJM.R turns the check back on where it tests it.
options(BJM.check_complete_case = FALSE)
