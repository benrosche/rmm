data <- 
  data.frame(pid=c(10,1,2,3,4,5,6,7,8,9,10,1,9,4,6,8,3,9,9,2,4), 
             gid=rep(1:10, times=c(2, 3, 1, 1, 2, 2, 2, 3, 2, 3))) %>%
  mutate(survival=rweibull(21, 1, 1),
         event=rbinom(21, 1, 0.5),
         var1=rnorm(21, 0, 1),
         var2=rnorm(21, 0, 1),
         var3=rnorm(21, 0, 1),
         var4=rnorm(21, 0, 1),
         var5=rnorm(21, 0, 1)) %>%
  group_by(gid) %>% add_count(name="var6")
