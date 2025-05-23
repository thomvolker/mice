% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{matchindex}
\alias{matchindex}
\title{Find index of matched donor units}
\usage{
matchindex(d, t, k = 5L)
}
\arguments{
\item{d}{Numeric vector with values from donor cases.}

\item{t}{Numeric vector with values from target cases.}

\item{k}{Integer, number of unique donors from which a random draw is made.
For \code{k = 1} the function returns the index in \code{d} corresponding
to the closest unit. For multiple imputation, the
advice is to set values in the range of \code{k = 5} to \code{k = 10}.}
}
\value{
An integer vector with \code{length(t)} elements. Each
element is an index in the array \code{d}.
}
\description{
Find index of matched donor units
}
\details{
For each element in \code{t}, the method finds the \code{k} nearest
neighbours in \code{d}, randomly draws one of these neighbours, and
returns its position in vector \code{d}.

Fast predictive mean matching algorithm in seven steps:
\enumerate{
\item Shuffle records to remove effects of ties
\item Obtain sorting order on shuffled data
\item Calculate index on input data and sort it
\item Pre-sample vector \code{h} with values between 1 and \code{k}
}

For each of the \code{n0} elements in \code{t}:
\enumerate{
\item find the two adjacent neighbours
\item find the \code{h_i}'th nearest neighbour
\item store the index of that neighbour
}

Return vector of \code{n0} positions in \code{d}.

We may use the function to perform predictive mean matching under a given
predictive model. To do so, specify both \code{d} and \code{t} as
predictions from the same model. Suppose that \code{y} contains the observed
outcomes of the donor cases (in the same sequence as \code{d}), then
\code{y[matchindex(d, t)]} returns one matched outcome for every
target case.

See \url{https://github.com/amices/mice/issues/236}.
This function is a replacement for the \code{matcher()} function that has
been in default in \code{mice} since version \code{2.22} (June 2014).
}
\examples{
set.seed(1)

# Inputs need not be sorted
d <- c(-5, 5, 0, 10, 12)
t <- c(-6, -4, 0, 2, 4, -2, 6)

# Index (in vector a) of closest match
idx <- matchindex(d, t, 1)
idx

# To check: show values of closest match

# Random draw among indices of the 5 closest predictors
matchindex(d, t)

# An example
train <- mtcars[1:20, ]
test <- mtcars[21:32, ]
fit <- lm(mpg ~ disp + cyl, data = train)
d <- fitted.values(fit)
t <- predict(fit, newdata = test)  # note: not using mpg
idx <- matchindex(d, t)

# Borrow values from train to produce 12 synthetic values for mpg in test.
# Synthetic values are plausible values that could have been observed if
# they had been measured.
train$mpg[idx]

# Exercise: Create a distribution of 1000 plausible values for each of the
# twelve mpg entries in test, and count how many times the true value
# (which we know here) is located within the inter-quartile range of each
# distribution. Is your count anywhere close to 500? Why? Why not?
}
\author{
Stef van Buuren, Nasinski Maciej, Alexander Robitzsch
}
