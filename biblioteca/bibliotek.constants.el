;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;                                                          ;;;;;;;;
;;;;;;  All files in this directory are copyright 1997, 1998    ;;;;;;;;
;;;;;;  by Rafael D. Sorkin.     All rights reserved.           ;;;;;;;;
;;;;;;                                                          ;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

        ;;; Definitions of some constants ;;;

;:----------------------------------------------
;:   LOAD THIS FILE BUT DO NOT COMPILE IT   
;:----------------------------------------------

;:                        ; Time-stamp:<21:49:23 May 30 1998 (13680 46883)>

;: These constant definitions are isolated here to avoid compilation, 
;  which causes trouble for both elisp and sun-lisp, at least. 
;
;  elisp trouble: the use of defconstant as an alias for defconst 
;  causes some weird error when compiling, but not when loading. 
;  (This will be fixed when elisp compiler learns how to alias a macro)
;
;  sun-lisp trouble: the use of fractions like 2/3 as variable names
;  discombobulates the sun-lisp compiler, because it "sees thru"
;  the `(when *elisp* ..)' wrapping. 

;===================================================================
(in-package 'user)                      ; this to make SCL happy
;===================================================================

;: Define some numerical constants

(defconstant the-number-pi 
  3.1415926535897932384626433832795028841972      " Pi to 40 places")
(defconstant pi%  the-number-pi                   " Pi to 40 places")

(defconstant the-number-e 
  2.7182818284590452353602874713526624977572      " e to 40 places")
(defconstant e%  the-number-e                     " e to 40 places")

(defconstant Euler% 
  0.5772156649015328606065120900824024310422  " Euler's constant to 40 places")

(defconstant log_sqrt_2pi  (/ (log (* 2 pi%)) 2))

;: Define some small fractions for elisp (which lacks rationals)
(when *elisp* 
  (setq  
       1/2 (/  1  2.0)                 
      -1/2 (/ -1  2.0)
       3/2 (/  3  2.0)
      -3/2 (/ -3  2.0)
       1/3 (/  1  3.0)
      -1/3 (/ -1  3.0)
       2/3 (/  2  3.0)
      -2/3 (/ -2  3.0)
       4/3 (/  4  3.0)
       1/4 (/  1  4.0)
       3/4 (/  3  4.0)
       1/5 (/  1  5.0)
       1/6 (/  1  6.0)       
       1/7 (/  1  7.0)
       1/8 (/  1  8.0)
       1/9 (/  1  9.0)
      1/10 (/  1 10.0)
      1/11 (/  1 11.0)
      1/12 (/  1 12.0)
    -1/360 (/ -1 360.0)
    1/1260 (/  1 1260.0)
   -1/1680 (/ -1 1680.0)
   ))
