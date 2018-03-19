%%%%%%%%%%%%%%%%%
%%%  Summary  %%%
%%%%%%%%%%%%%%%%%

The goal of the codes in this project is to develop a numerical implementation of quantum corrections in plasmonics. This work is based on the ideas in "Christensen, T., Yan, W., Jauho, A. P., Soljačić, M., & Mortensen, N. A. (2017). Quantum corrections in nanoplasmonics: shape, scale, and material. Physical Review Letters, 118(15), 157402". The factors that are used with the Feibelman d-parameters to compute quantum corrections are called Lambdas. These Lambdas are the objective of these functions.

These codes are intended to be used along with the MNPBEM Toolbox, by Ulrich Hohenester (http://physik.uni-graz.at/mnpbem/), since we make use some of their functions and solvers.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Compatibility with MNPBEM	%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

These functions cannot be used on their own, they need the MNPBEM directory to be in the Matlab path.

IMPORTANT: to successfully use the colormap 'bluewhitered_mod()', which sets white for 0 and red (blue) for positive (negative) values, when displaying the surface charge of every mode of a system: please add the following line at the end of the file MNPBEM14>Misc>@bemplot>refresh.m:
>>colormap(bluewhitered_mod()) 

%%%%%%%%%%%%%%
%%%  Help  %%%
%%%%%%%%%%%%%%

Each function contains an explanation, as well as information about the inputs and outputs.
Some examples on how to use the functions can be found in the examples folder.
              
%%%%%%%%%%%%%%%%%%%
%%%  Copyright  %%%
%%%%%%%%%%%%%%%%%%%
MIT License

Copyright (c) 2018 Álvaro Gómez Iñesta.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

 