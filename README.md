# lmld_amnat_2015
Code for spatial population dynamics from Usinowicz 2015
%Notes Dec5
%Using this code mostly to mess around with approximations at the very end.
%They all give close fits. Two versions of the linearized, low-density
%version: simple quadratics (Ni=4*x+4 or 2*x+6, Nm=6*x+2), full
%approximation (Ni=a*x+p, Nm=b*x+p). There may still be some questions with
%Nr: does it work best with x, c*x+p-p. Is c based on 1 or 1/2? By eye,
%best fit right now is with full approximations, Nr=x. 


% This code is a cleaned up version of the most current code tht I am
% using to run the spatial lottery model -- so no stages. Look at the code
% in usinowicz_3stg_FDM.m for a non-spatial, stage-structure model. I've
% tried to make notes throughout to explain what things are. 

%The basic structure of this program is three large iterative loops. The
%two outermost loops have nothing to do with the mechanics of the model,
%they are simply for doing reps to build some stats, or to increment
%through certain parameters (in this program, the size of the starting
%cluster). You can basically ignore them for learning how the model works. 

%There are actually two model imbedded in this code -- the spatial (variable pop)
%and non-spatial (variable pop_lot) versions of the lottery model. The 
%non-spatial (classic) model is really simple and could probably be written in 
%about 12 lines of code. I've tried to comment around it to show you where it is. 

%Most of the code is for the spatial model. The most difficult aspect of
%this code to understand is probably the dispersal step. The easiest way to
%do this in Matlab is to make use of its built-in Fourier transform
%algorithms and to do a convultion of the matrix of individuals, with a
%matrix defining the dispersal kernel (in Fourier space). This is a very
%mathematically complicated thing, but a very easy tool to use in
%Matlab. If you want to learn more about it, I can dig up some web sites
%that I found initially that did a good job explaining this stuff. 
