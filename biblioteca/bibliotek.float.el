;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;                                                          ;;;;;;;;
;;;;;;  All files in this directory are copyright 1997, 1998    ;;;;;;;;
;;;;;;  by Rafael D. Sorkin.     All rights reserved.           ;;;;;;;;
;;;;;;                                                          ;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


                            ; Time-stamp:<23:52:41 21-Jun-1997 (13228 44953)>

; Functions involved with floating point computations, especially ones where
; the result produced is inexact for reasons beyond truncation error of the
; result itself. 
; 
;
; Numerical constants defined 
; ===========================
;
;  These constants are all defined in a special file because they give trouble
;  when compiling!   They are:  
;
;    pi, e, Euler, (log (sqrt 2 pi)),  
;    and a bunch of fractions like 1/3 for elisp
;
;  Infinity NaN are defined in preparations.el
;
; Functions defined (or waiting to be defined)
; ============================================
;
;  sinh cosh tanh
;
;  polynomial-value
;
;  !                   (factorial of integer argument -- elisp version)
;                       (It's in this bibliotek because it converts to float,
;                        TCL version is in bibliotek.general)
;  log-factorial       (fairly accurate, useful for large arguments)
;  factorial           (for now just relies on log-factorial)
;
;  solve-increasing-fcn  (solve y=f(x) for x) (presently in developing/)
;
;  A bunch of matrix functions
;
;    For now these are in a separate file
;
;=======================================================================
(in-package 'user)                      ; merely to make SCL happy
;=======================================================================

(when *elisp*                           ; SCL has them already
  (defun sinh(x) 
    (let ((y (exp x)))
      (/ (- y (reciprocal y)) 2)))
  (defun cosh(x) 
    (let ((y (exp x)))
      (/ (+ y (reciprocal y)) 2)))
  (defun tanh(x)
    (let* ((y (exp x))
           (z (reciprocal y)))
      (/ (- y z) (+ y z)))))

;------------------------------------------------------------------------

(defun polynomial-value (p x)
  "\
 The first arg should be the LIST of coefficients  (a_0 a_1...a_n), 
 with 0's for absent terms, of course.
 NB: the coefficients are in order of INCREASING DEGREE"
  (if (null p) 0
    (+
     (car p)
     (* x (polynomial-value (cdr p) x)))))
 ;
 ; This is the more accurate and efficient way to evaluate a polynomial.
 ; Seems no need to localize the symbols.

;--------------------------------------------------------------------------
 
(deff log-factorial (x &optional (thresh 20.0))  
  "\
 Natural log of factorial of x, where x is a real number > -1.
 Evaluation is done using Stirling series, after first promoting x to
 be greater than the threshold specified by the optional second argument
 (which defaults to 20).
  "
  (cond
  ;-----------------------
  ; check that arg is > -1
  ;-----------------------
   ((<= x -1)
    (error "`log-factorial' requires arg > -1 for now."))
  ;------------------------------------------------------------------
  ; handle x=0 and 1 specially to avoid roundoff error in those cases
  ;------------------------------------------------------------------
   ((or 
     (= x 0)
     (= x 1))
              0.0)
  ;------------------------------------------
  ; if x < thresh then augment it and recurse
  ;------------------------------------------
   ((< x thresh)
    (- (log-factorial (1+ x) thresh) (o log 1+ x)))  
  ;------------------------------------
  ; if x >= thresh then plug into the series
  ;------------------------------------
   (t (+
       (* (+ x 1/2) (log x))
       (- x)
       log_sqrt_2pi
       (polynomial-value 
         (list 0 1/12 0 -1/360 0 1/1260 0 -1/1680)
         (reciprocal x))))))
         ;
  ; The method is to make the argument large enough that a few terms in the
  ; Stirling series will be quite accurate. 
  ;
  ; It seems to work well, but behaves rather strangely.
  ; First of all, one might have expected the optimum threshold to be
  ; around 33 (because there the ratio of the smallest to biggest term is
  ; about 10^16) but empirically 20 appears to work best (see the tests).
  ; Moreover, many choices of threshshold in the range 15 to 80 (if not
  ; beyond) give precisely the same answer (see the tests).
  ;
  ; Improvement
  ;  Let it handle args < -1 as well.  One nuisance though, is that log is
  ;  complex in parts of this range, so maybe should give just its real part
  ;  (i.e. (log abs fact x))
  ;
  ; The coefficients in the series come from Dwight 851.5 (page 210)
  ;  Here is how they were computed
  ;
  ; (* 6  1 2) ; => 12
  ; (* 30 3 4) ; => 360
  ; (* 42 5 6) ; => 1260
  ; (* 30 8 7) ; => 1680
  ;
  ; (- (log (sqrt (atl 2 * pi))) log_sqrt_2pi) ; => 0.0

;; Here is older version, with thresh built in as 20, rather than being an
;; optional parameter
;
; (defun log-factorial (x)  
;   " Natural log of x! using Stirling series."
;   (cond
;    ((= x 0) 0.0)
;    ((< x 20) (- (log-factorial (1+ x)) (log (1+ x))))  ; Notice choice of 20
;    (t (+
;	(* (+ x 1/2) (log x))
;	(- x)
;	log_sqrt_2pi
;	(polynomial-value 
;	  (list 0 1/12 0 -1/360 0 1/1260 0 -1/1680)
;	  (reciprocal x))))))

;--------------------------------------------------------------------------

(defun factorial (x) 
  "\
 A poor excuse for this function scrabbled together from wherever.
 Among other problems, it rejects arguments < -1"
  (if (and (integerp x) (< x 171)) (! x)
    (o exp log-factorial x)))

;--------------------------------------------------------------------------

(when *elisp* 
(defun ! (n) 
  "\
 Factorial function for integers -- from 0 to 170 only!  
 Converts argument to float to avoid large integer problem in elisp."
  (cond 
   ((not (integerp n)) (error "This one is for integers only"))
   ((< n 0) Infinity%)
   ((= 0 n) 1.0)
   (t (* n (! (1- n))))))  )
  ;
  ; we make the result float to avoid the large integer problem in elisp
  ; The TCL version is better and is in bibliotek.general

;--------------------------------------------------------------------------
