app <- ShinyDriver$new("../../")
app$snapshotInit("ex1")

app$setInputs(lambda_r = "0.677, 0.746, 0.768, 0.664")
app$setInputs(tau_r = "2.275, 2.928, 3.437, 3.095")
app$setInputs(usetau_f = TRUE)
app$setInputs(tau_f = "2.275, 2.928, 3.013, 3.095")
app$setInputs(theta_r = "0.744, 1.061, 0.899, 0.638")
app$setInputs(kappa_r = -0.032)
app$setInputs(usekappa_f = TRUE)
app$setInputs(kappa_f = 0)
app$setInputs(phi_r = 1.032)
app$setInputs(usephi_f = TRUE)
app$setInputs(phi_f = 1)
app$setInputs(usepropsel = TRUE)
app$setInputs(prop = 0.25)
app$snapshot()
