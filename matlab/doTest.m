
clear;
clc;

%site_points = [6 7; 7 6];

%site_points = [ 7 7; 5 6; 6 2 ];

%site_points = [ 7 7; 5 6; 6 2; 8 3];

site_points = [ 7 7; 4 12; 5 10; 9 14; 10 9; 3 9; 6 20; 15 15; 12 13];

v = VoronoiFortuneAlgo(site_points, 0.5);

v.do();
