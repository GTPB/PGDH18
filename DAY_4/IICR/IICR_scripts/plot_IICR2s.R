n <- 10
M <- 1
t <- seq(0.001, 5, 0.001)
IICR_t <- IICR2nislands(n, M, t)
plot(t, IICR_t)