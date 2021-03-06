% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loxcode_sample.R
\docType{methods}
\name{impute,loxcode_sample-method}
\alias{impute,loxcode_sample-method}
\title{Impute missing code in 13-element cassettes}
\usage{
\S4method{impute}{loxcode_sample}(x)
}
\description{
For 13-element cassettes that are missing a single element, the
missing element is imputed to minimise the resulting dist_orig.
}
