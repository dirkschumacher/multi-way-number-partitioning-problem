
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Multi-way-number-partitioning

<!-- badges: start -->

<!-- badges: end -->

## Introduction

The multi-way number partitioning problem tries to partition the
elements of a finite set \(I\) of integers into \(K > 0\) mutually
exclusive subsets \(I_k\) such that their sum is as close together as
possible. More formally: let \(S_k = \sum_{i \in I_k}i\), then the goal
is to minimize \(\max_{k \in K}S_k - \min_{k \in K}S_k\)
[ref](https://en.wikipedia.org/wiki/Partition_problem).

## Mixed-integer programming model

One way to solve such a problem is to formulate it as a mixed-integer
linear program.

### Variables

For each element \(i\) and \(k\) \(x_{i,k}\) is 1 if element \(i\) is
part of set \(S_k\).

``` r
x <- model$add_variable(
  i = 1:n_I,
  k = 1:K,
  type = "integer",
  lb = 0,
  ub = 1
)
```

Then we need three more variables to keep track of the size of the
subsets.

``` r
size <- model$add_variable(k = 1:K)
max_size <- model$add_variable()
min_size <- model$add_variable()
```

### Objective

Minimize the difference between the largest and the smallest set.

``` r
model$set_objective(max_size - min_size, sense = "min")
```

### Constraints

Each element \(i\) needs to part of exactly one subset \(S_k\).

``` r
model$add_constraint(
  sum_expr(x[i, k], k = 1:K) == 1,
  i = 1:n_I
)
```

The following constraints model the correct size of each subset and
define the correct values of the minimum and maximum subset values.

``` r
model$add_constraint(
  sum_expr(x[i, k] * I[i], i = 1:n_I) == size[k],
  k = 1:K
)
model$add_constraint(size[k] >= min_size, k = 1:K)
model$add_constraint(size[k] <= max_size, k = 1:K)
```

Some notes:

  - I did test a formulation which orders the subsets by size to break
    symmetries. This leads to worse computational results compared to
    the model above. Maybe for larger numbers of \(K\) this makes sense.
  - However I do need to run a full computational study with a larger
    number of samples to make any claims.
  - Once I have a larger testset I will also test CBC.
  - This needs to be compared with the greedy heuristic implemented by
    Noam Ross.

### Computational results

We use *GLPK* as a solver and see how far we can go.

#### The implementation

``` r
library(rmpk)
library(ROI.plugin.glpk)
solve_partitioning_problem <- function(I, K, time_limit_sec = 30) {
  n_I <- length(I)
  model <- MIPModel(ROI_solver("glpk", control = list(verbose = TRUE, presolve = TRUE, tm_limit = time_limit_sec * 1000)))
  x <- model$add_variable(
    i = 1:n_I,
    k = 1:K,
    type = "integer",
    lb = 0,
    ub = 1
  )
  size <- model$add_variable(k = 1:K)
  max_size <- model$add_variable()
  min_size <- model$add_variable()
  model$set_objective(max_size - min_size, sense = "min")
  model$add_constraint(
    sum_expr(x[i, k], k = 1:K) == 1,
    i = 1:n_I
  )
  model$add_constraint(
    sum_expr(x[i, k] * I[i], i = 1:n_I) == size[k],
    k = 1:K
  )
  model$add_constraint(size[k] >= min_size, k = 1:K)
  model$add_constraint(size[k] <= max_size, k = 1:K)
  model$optimize()
  list(model = model, solution = model$get_variable_value(x[i, k]))
}
```

We also need a function score a solution and to display it.

``` r
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))
display <- function(result, I) {
  solution <- result$solution
  solution %>% 
    filter(value == 1) %>% 
    group_by(k) %>% 
    summarise(sum = sum(I[i]), n = n())
}
score <- function(result, I) {
  subset_sum <- display(result, I)$sum
  max(subset_sum) - min(subset_sum)
}
```

#### Small example

``` r
set.seed(42)
K <- 4
n_I <- 10
I <- sample.int(n_I, size = 10000, replace = TRUE, prob = runif(n_I))
I <- sort(as.integer(table(I)), decreasing = TRUE)
```

``` r
result <- solve_partitioning_problem(I, K)
#> <SOLVER MSG>  ----
#> GLPK Simplex Optimizer, v4.65
#> 22 rows, 46 columns, 100 non-zeros
#> Preprocessing...
#> 22 rows, 46 columns, 100 non-zeros
#> Scaling...
#>  A: min|aij| =  1.000e+00  max|aij| =  1.550e+03  ratio =  1.550e+03
#> GM: min|aij| =  9.301e-01  max|aij| =  1.075e+00  ratio =  1.156e+00
#> EQ: min|aij| =  8.894e-01  max|aij| =  1.000e+00  ratio =  1.124e+00
#> Constructing initial basis...
#> Size of triangular part is 22
#>       0: obj =   0.000000000e+00 inf =   5.469e+02 (1)
#>       5: obj =   5.000000000e+03 inf =   1.776e-14 (0)
#> *    16: obj =  -4.547473509e-13 inf =   7.461e-14 (0)
#> OPTIMAL LP SOLUTION FOUND
#> GLPK Integer Optimizer, v4.65
#> 22 rows, 46 columns, 100 non-zeros
#> 40 integer variables, all of which are binary
#> Preprocessing...
#> 22 rows, 46 columns, 100 non-zeros
#> 40 integer variables, all of which are binary
#> Scaling...
#>  A: min|aij| =  1.000e+00  max|aij| =  1.550e+03  ratio =  1.550e+03
#> GM: min|aij| =  9.301e-01  max|aij| =  1.075e+00  ratio =  1.156e+00
#> EQ: min|aij| =  8.894e-01  max|aij| =  1.000e+00  ratio =  1.124e+00
#> 2N: min|aij| =  5.000e-01  max|aij| =  1.625e+00  ratio =  3.250e+00
#> Constructing initial basis...
#> Size of triangular part is 22
#> Solving LP relaxation...
#> GLPK Simplex Optimizer, v4.65
#> 22 rows, 46 columns, 100 non-zeros
#>      16: obj =  -1.000000000e+04 inf =   2.500e+03 (4)
#>      35: obj =  -4.547473509e-13 inf =   2.961e-14 (0)
#> OPTIMAL LP SOLUTION FOUND
#> Integer optimization begins...
#> Long-step dual simplex will be used
#> +    35: mip =     not found yet >=              -inf        (1; 0)
#> +    67: >>>>>   3.330000000e+02 >=   0.000000000e+00 100.0% (17; 0)
#> +   203: >>>>>   2.810000000e+02 >=   0.000000000e+00 100.0% (51; 12)
#> +   225: >>>>>   1.390000000e+02 >=   0.000000000e+00 100.0% (56; 24)
#> + 11647: mip =   1.390000000e+02 >=     tree is empty   0.0% (0; 7253)
#> INTEGER OPTIMAL SOLUTION FOUND
#> <!SOLVER MSG> ----
```

``` r
stopifnot(result$model$termination_status() == rmpk::TERMINATION_STATUS$SUCCESS)
```

Let’s check the score. Note: this should be the same value as the final
upper bound in the log file.

``` r
score(result, I)
#> [1] 139
```

``` r
stopifnot(score(result, I) == result$model$objective_value())
```

``` r
display(result, I)
#> # A tibble: 4 x 3
#>       k   sum     n
#>   <int> <int> <int>
#> 1     1  2507     2
#> 2     2  2418     2
#> 3     3  2557     3
#> 4     4  2518     3
```

#### Small K, large I

``` r
set.seed(42)
K <- 4
n_I <- 200
I <- sample.int(n_I, size = 3e6, replace = TRUE, prob = runif(n_I))
I <- sort(as.integer(table(I)), decreasing = TRUE)
```

``` r
result <- solve_partitioning_problem(I, K)
#> <SOLVER MSG>  ----
#> GLPK Simplex Optimizer, v4.65
#> 212 rows, 806 columns, 1620 non-zeros
#> Preprocessing...
#> 212 rows, 806 columns, 1620 non-zeros
#> Scaling...
#>  A: min|aij| =  1.000e+00  max|aij| =  2.871e+04  ratio =  2.871e+04
#> GM: min|aij| =  9.359e-01  max|aij| =  1.068e+00  ratio =  1.142e+00
#> EQ: min|aij| =  8.984e-01  max|aij| =  1.000e+00  ratio =  1.113e+00
#> Constructing initial basis...
#> Size of triangular part is 212
#>       0: obj =   0.000000000e+00 inf =   5.090e+04 (1)
#>      63: obj =   1.500000000e+06 inf =   3.341e-13 (0)
#> *   196: obj =   3.492459655e-10 inf =   4.833e-13 (0)
#> OPTIMAL LP SOLUTION FOUND
#> GLPK Integer Optimizer, v4.65
#> 212 rows, 806 columns, 1620 non-zeros
#> 800 integer variables, all of which are binary
#> Preprocessing...
#> 212 rows, 806 columns, 1620 non-zeros
#> 800 integer variables, all of which are binary
#> Scaling...
#>  A: min|aij| =  1.000e+00  max|aij| =  2.871e+04  ratio =  2.871e+04
#> GM: min|aij| =  9.359e-01  max|aij| =  1.068e+00  ratio =  1.142e+00
#> EQ: min|aij| =  8.984e-01  max|aij| =  1.000e+00  ratio =  1.113e+00
#> 2N: min|aij| =  5.000e-01  max|aij| =  1.366e+00  ratio =  2.732e+00
#> Constructing initial basis...
#> Size of triangular part is 212
#> Solving LP relaxation...
#> GLPK Simplex Optimizer, v4.65
#> 212 rows, 806 columns, 1620 non-zeros
#>     196: obj =  -3.000000000e+06 inf =   1.875e+05 (4)
#>     319: obj =   1.500000000e+06 inf =   0.000e+00 (0)
#> *   524: obj =   0.000000000e+00 inf =   3.013e-12 (0)
#> OPTIMAL LP SOLUTION FOUND
#> Integer optimization begins...
#> Long-step dual simplex will be used
#> +   524: mip =     not found yet >=              -inf        (1; 0)
#> +   997: >>>>>   1.356000000e+03 >=   0.000000000e+00 100.0% (294; 0)
#> +  2266: >>>>>   9.740000000e+02 >=   0.000000000e+00 100.0% (939; 19)
#> +  2528: >>>>>   6.880000000e+02 >=   0.000000000e+00 100.0% (1076; 68)
#> +  3832: >>>>>   4.580000000e+02 >=   0.000000000e+00 100.0% (1361; 414)
#> +  5513: >>>>>   3.880000000e+02 >=   0.000000000e+00 100.0% (1621; 795)
#> +  6007: >>>>>   3.470000000e+02 >=   0.000000000e+00 100.0% (1755; 867)
#> + 11107: >>>>>   2.500000000e+02 >=   0.000000000e+00 100.0% (3695; 988)
#> + 12449: >>>>>   2.030000000e+02 >=   0.000000000e+00 100.0% (3105; 3208)
#> + 18408: >>>>>   1.580000000e+02 >=   0.000000000e+00 100.0% (4886; 3349)
#> + 21449: >>>>>   1.450000000e+02 >=   0.000000000e+00 100.0% (6005; 3428)
#> + 21537: >>>>>   1.430000000e+02 >=   0.000000000e+00 100.0% (6029; 3434)
#> + 24736: >>>>>   1.230000000e+02 >=   0.000000000e+00 100.0% (6991; 3492)
#> + 36834: mip =   1.230000000e+02 >=   0.000000000e+00 100.0% (9465; 5807)
#> + 39650: mip =   1.230000000e+02 >=   0.000000000e+00 100.0% (10360; 5848)
#> + 45469: mip =   1.230000000e+02 >=   0.000000000e+00 100.0% (11939; 5933)
#> + 53560: mip =   1.230000000e+02 >=   0.000000000e+00 100.0% (13921; 6034)
#> TIME LIMIT EXCEEDED; SEARCH TERMINATED
#> <!SOLVER MSG> ----
```

After 30 seconds the score is: 123

``` r
display(result, I)
#> # A tibble: 4 x 3
#>       k    sum     n
#>   <int>  <int> <int>
#> 1     1 750088    47
#> 2     2 749975    64
#> 3     3 749965    41
#> 4     4 749972    48
```
