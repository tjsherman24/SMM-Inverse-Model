TransitionMatrixCalculator
%TransitionMatrixCalculator takes time data from Breakthrough Curves 1&2
%(tau1,tau2) as well as correpsonding concentration data (count1,count2)
%and estimates the transition Matrix between curves.  

%This script completes the 'Finding PDFs of arrival times', 'Classifying
%Arrival Times' and 'Construction of Systems of Equations' in our
%correspoding paper

%Inputs
%tau1,tau2: time data for BTC1, BTC2
%count1,count2: concentration data for BTC1,BTC2

%User Variables
%bn - set number of transition classes
%dec - sets how to round arrival times. dec 3 discretizing times to nearest 1000th place-see script for detail


UncertaintyMatrices
%UncertaintyMatrices estimates the mean, lower bound, and upper bound
%transition matrices.  The user can specify what values to filter out, i.e.
%if estimated matrix elements have unphysical values greater than 1, you
%can apply a filter. Default filter, rids elements greater than .5

% The matrix is divided into different sections. Each corner  is a 4x4
% block and corresponds to a section. We also have sections that correspond
% to the middle rows and middle columns. For each section, 'b' corresponds to
% minimum filter, meaning all estimated values below this value will be
% filtered out and 'c' correponds to the max filter, meaning all values
% above this threshold will be filtered. 

% This script Solves the systems of equations and calculates Upper and
% Lower bounded matrices, as well as mean matrix

%Outputs:
%M_guess - mean matrix
%M_high - upper bound matrix
%M_low - lower bound matrix