data("coalgov")

test_that("dissectFormula() works with all model families and with and without intercept", {
  
  # Gaussian
  expect_no_error(
    dissectFormula(
      formula(sim.y ~ 1 + majority + mwc + hm(id = cid, name = cname, type = RE, showFE = F) + mm(id(pid, gid), mmc(ipd), mmw(w ~ 1/offset(n), constraint = 1))),
      family="Gaussian",
      data
    ),
    message="LHS CHECK (Gaussian)"
  )
  
  expect_error(
    dissectFormula(
      formula(Surv(sim.st, sim.e) ~ 1 + majority + mwc + hm(id = cid, name = cname, type = RE, showFE = F) + mm(id(pid, gid), mmc(ipd), mmw(w ~ 1/offset(n), constraint = 1))),
      family="Gaussian",
      data
    ),
    label="LHS CHECK (Gaussian)"
  )
  
  
  # Binomial 
  expect_no_error(
    dissectFormula(
      formula(sim.y ~ 1 + majority + mwc + hm(id = cid, name = cname, type = RE, showFE = F) + mm(id(pid, gid), mmc(ipd), mmw(w ~ 1/offset(n), constraint = 1))),
      family="Binomial",
      data
    ),
    message="LHS CHECK (Binomial)"
  )
  
  expect_error(
    dissectFormula(
      formula(Surv(sim.st, sim.e) ~ 1 + majority + mwc + hm(id = cid, name = cname, type = RE, showFE = F) + mm(id(pid, gid), mmc(ipd), mmw(w ~ 1/offset(n), constraint = 1))),
      family="Binomial",
      data
    ),
    label="LHS CHECK (Binomial)"
  )

  # Weibull
  expect_no_error(
    dissectFormula(
      formula(Surv(sim.st, sim.e) ~ 1 + majority + mwc + hm(id = cid, name = cname, type = RE, showFE = F) + mm(id(pid, gid), mmc(ipd), mmw(w ~ 1/offset(n), constraint = 1))),
      family="Weibull",
      data
    ),
    message="LHS CHECK (Weibull)"
  )
  
  expect_error(
    dissectFormula(
      formula(sim.y ~ 1 + majority + mwc + hm(id = cid, name = cname, type = RE, showFE = F) + mm(id(pid, gid), mmc(ipd), mmw(w ~ 1/offset(n), constraint = 1))),
      family="Weibull",
      data
    ),
    label="LHS CHECK (Weibull)"
  )
  
  expect_no_error(
    dissectFormula(
      formula(sim.y ~ majority + mwc + hm(id = cid, name = cname, type = RE, showFE = F) + mm(id(pid, gid), mmc(ipd), mmw(w ~ 1/offset(n), constraint = 1))),
      family="Gaussian",
      data
    ),
    message="INTERCEPT MISSING"
  )
  
  # 2do: add Cox
  
})

test_that("dissectFormula() works with different ways of specifying mm()", {
  
  # RHS

  expect_error(
    dissectFormula(
      formula(Surv(govdur, earlyterm) ~ 1 + mm(mmc(ipd), mmw(w ~ 1/exp(ipd), c = 2))),
      family="Weibull",
      data
    ),
    label="mm() -> id() missing"
  )
  
  expect_error(
    dissectFormula(
      formula(Surv(govdur, earlyterm) ~ 1 + mm(id(pid, gid), mmw(w ~ 1/exp(ipd), c = 2))),
      family="Weibull",
      data
    ),
    label="mm() -> mmc() missing"
  )
  
  expect_error(
    dissectFormula(
      formula(Surv(govdur, earlyterm) ~ 1 + mm(id(pid, gid), mmc(ipd))),
      family="Weibull",
      data
    ),
    label="mm() -> mmw() missing"
  )
  
  expect_no_error(
    dissectFormula(
      formula(Surv(govdur, earlyterm) ~ 1 + mm(id(pid, gid), mmc(), mmw(w ~ 1/offset(n), constraint = 1))),
      family="Weibull",
      data
    ),
    message="mm() -> no mmc() variables should be possible"
  )
  
  expect_equal(
    {dissectFormula(
      formula(Surv(govdur, earlyterm) ~ 1 + mm(id(pid, gid), mmc(), mmw(w ~ 1/offset(n), constraint = 1))),
      family="Weibull",
      data)}[["l1"]][["mmwconstraint"]], 
    1, 
    label="mm() -> constraint = 1 should return 1"
  )
  
  expect_equal(
    {dissectFormula(
      formula(Surv(govdur, earlyterm) ~ 1 + mm(id(pid, gid), mmc(), mmw(w ~ 1/offset(n), c = 2))),
      family="Weibull",
      data)}[["l1"]][["mmwconstraint"]], 
    2, 
    label="mm() -> c = 2 should return 2"
  )
  
  expect_equal(
    {dissectFormula(
      formula(Surv(govdur, earlyterm) ~ 1 + mm(id(pid, gid), mmc(), mmw(w ~ 1/offset(n)))),
      family="Weibull",
      data)}[["l1"]][["mmwconstraint"]], 
    1, 
    label="mm() -> no constraint should return 1"
  )
  
})

test_that("dissectFormula() works with different ways of specifying hm()", {
  
  expect_error(
    dissectFormula(
      formula(sim.y ~ 1 + hm(name = cname, type = RE, showFE = F)),
      family="Gaussian",
      data
    ),
    label="hm() -> id() missing"
  )
  
  expect_error(
    dissectFormula(
      formula(Surv(govdur, earlyterm) ~ 1 + hm(id = cid, name =, type = RE, showFE = F)),
      family="Weibull",
      data
    ),
    label="hm() -> name improperly specified"
  )
  
  expect_no_error(
    dissectFormula(
      formula(Surv(govdur, earlyterm) ~ 1 + hm(id = cid, name = cname)),
      family="Weibull",
      data
    ),
    message="hm() -> no type or showFE should be possible"
  )
  
 
  
})

