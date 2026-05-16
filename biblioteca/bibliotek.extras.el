;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;                                                          ;;;;;;;;
;;;;;;  All files in this directory are copyright 1997, 1998,   ;;;;;;;;
;;;;;;  1999 by Rafael D. Sorkin.  All rights reserved.         ;;;;;;;;
;;;;;;                                                          ;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;                              Time-stamp:<Nov 03 1998 02:46:57 (13886 46321)>

;: This is the file { bibliotek.extras.el }
;: Its contents are here for temporary convenience,
;: intended to go into other biblioteks eventually, as indicated.


;: Roster of functions/macros

 ;  make-KR-preorder (--> poset)
 ;
 ;  make-test-preposet (--> poset)
 ;
 ;  isomorphism-p (--> poset)
 ;
 ;  compiled-eval (evals after compiling) (--> macros)
 ;
 ;  choose  (binomial coefficient)  (--> general)
 ;
 ;  sup-norm (--> general)
 ;
 ;  menge-p/setp (--> general)
 ; 
 ;  bijection-p (--> general)
 ;
 ;  solve-monotone (--> general)
 ;
 ;  list-plist    (nicely prints the entries in the plist L) (--> general)
 ;
 ;  print-to-file  (--> general)
 ;  print-list-as-line  (prints multiple args on a single line)
 ;  print-list-in-lines (prints multiple args on multiple lines)
 ;  print-as-line       (prints single arg on its own line)
 ;  (print-plist) (alias for same)

;=================================================================
(in-package 'user)                      ; merely to make SCL happy
;=================================================================

;: The functions and macros

;:: Various ones 

(deff make-KR-preorder (N name)
  "\
 Make a 3 layer poset of the KR type with `N' (anonymous) elts and `name'
 as symbolic name. (Note: we do not impose that every elt of the bottom layer
 precedes every elt of the top; we just insert random links between adjacent
 layers) "
  (&localize N name R top-numb bottom-numb bottom middle top x y)
  (prepare-substrate :name name :N N)
  (varbind
   R (eval name)
   top-numb (Round N 4)
   bottom-numb (Round (- N top-numb) 3)
   bottom (subseq R 0 bottom-numb)
   middle (subseq R bottom-numb (- N top-numb))
   top (subseq R (- N top-numb) N))
  ;
  (loop for x in bottom do
  (loop for y in middle
    if (= (Random 2) 1)
    do (push x (past y))))
  ;
  (loop for x in middle do
  (loop for y in top
    if (= (Random 2) 1)
    do (push x (past y)))))


(deff make-test-preposet (&key 
			   ((:size N))
			   ((:p p_0))
			   ((:name name) 'P))
  "\
  The arguments are  
		   :size   :p   :name [default `P']
  We ``percolate'' 
  a preposet to use in tests of poset fcns.   The orders of the elts 
  in the substrate and in the individual pasts are all randomized.
  You *must* specify :size and :p.  
   "
  (&localize N p_0 name PS j k x)
  (varbind PS nil j nil k nil)
  (require 'preparations)
  (unless (and p_0 N)
    (error "Must specify both size and ``bare'' probability"))
  (prepare-disconnected-poset :name name :N N :anon nil)
  (setq PS (symbol-value name))
 ;-------------------------
 ; Generate random "bonds"
 ;-------------------------
  (kwa j from 0 upto N
  (kwa k from (1+ j) upto N
    (when (<  (Random 1.0) p_0) 
      (push 
       (nth j PS) 
       (past (nth k PS))))))
 ;-------------------------------------
 ; Randomly permute substrate and pasts
 ;-------------------------------------
  (set name (shuffle PS))
  (loop for x in PS do (setf (past x) (shuffle (past x))))
 ;-------------------------------------------
 ; Finally return the symbol naming the poset
 ;-------------------------------------------
  name)

(deff menge-p (S) 
  "\
 Is the argument a list bila repetitions? 
 Comparison done with `eq' NOT with `equal'.
  "
  (and
    (listp S)
    (loop 
     for cell on S
     if (memq (car cell) (cdr cell))
     return nil
     finally (return t))))
   ;
(defalias 'setp 'menge-p)


(deff bijection-p (f A B)
  "\
 Is f:A-->B a bijection?  This ONLY works if A and B are sets, 
 ie lists BILA duplications
  "
  (&localize f A B)
  (when *carefully*
    (unless (and (setp A)(setp B))
      (error "Argument to `bijection-p' is not a set")))
  (and
    (= (length A) (length B))
    (subsetp B (mapcar f A))))


(deff isomorphism-p (phi P Q)
  "\
 Are the orders P and Q isomorphic?  
 (If not put in *mrv* the first elts whose pasts fails to agree.)
 Actually P and Q can be any relations, they needn't be orders. \
  "
  (&localize phi P Q x)
  (makunbound '*mrv*)
  (and
   (bijection-p phi P Q)
   (loop 
     for x in P
     unless
       (equal-as-sets
	 (past (funcall phi x))
	 (mapcar phi (past x)))
     do 
       (setq *mrv* (list x (funcall phi x)))
       (message " elements x (phi x) have been placed in *mrv*")
     and return nil
     finally (return t))))


(LUN "compiled-eval" (H)
 (defmacro compiled-eval (&rest H)
  "evals the forms that follow after compiling them"
  `(funcall
    (Compile
     (lambda () ,@ H)))))


(deff choose (n &optional m &rest D)
  "\
 Generalized binomial coefficient for integer args (for now anyhow).
 The answer depends on the number of arguments:
  if 1 arg then return factorial
  if 2 args n m then return ``n choose m''
  if 3 or more args n a b ... c then return n!/(a! b!...c!)
 BEWARE: 
  In case of elisp answer may not be an exact integer if args are big.
  If it works for negative integer arguments, it's a bonus!
  "
  (cond
   ((not m) (! n))
   ;
   ((not D) 
    (varbind 
     n-m (- n m)
     a (max m n-m)
     b (min m n-m))
    (/
     (product (loop for k from (1+ a) to n collect k))
     (! b)))
   ;
   (t 
    (unless (= n (sum (cons m D)))
      (error "Invalid arguments to `choose'"))
    (*
     (choose n (- n m))
     (apply #'choose (cons (- n m) D))))))
 ;;
 ;; Should this just be called `!' itself, or should we reserve that to mean
 ;; n(n-1)(n-2)...(n-m+1)? 


(defun sup-norm (Y) 
  "\
 Maximum absolute value of an element of Y.  Here Y can be a single number,
 a sequence of numbers (eg a vector), or in general any sequence
 objects implied inductively by these.   
 (This includes a matrix, if it is implemented as a vector of vectors, but
 we'd have to extend if we wanted an array in TCL.)  \
  "
  (cond
   ((numberp Y) (abs Y))
   ((stringp Y) 
    (error "Character string given to `sup-norm'"))   ; only elisp has this
						      ; problem  
   ((sequencep Y) (sup (map 'list 'sup-norm Y)))
   (t 
    (error "Bad data to function `sup-norm'"))))

(deff solve-monotone (f y &key 
			   ((:ii ii) (list -1 1))
			   ((:tol-x tol-x) 1e-15) 
			   ((:tol-y tol-y) 1e-14)
			   ((:maxit maxit) 512))
  "\
 Solves the equation y=f(x) for x, where f is a monotonically 
 increasing or decreasing function.  Arguments are:
   1. the function f (or a symbol for it)
   2. the target value y
   3. :ii = (a b) is the `initial interval' in which to begin the search
            (defaults to [-1 1]) 
   4. :tol-x = (absolute) tolerance for x (defaults to 1e-15)
   5. :tol-y = (absolute) tolerance for y (defaults to 1e-14)
   6. :maxit = maximum number of iterations before quitting (defaults to 512)"
 ;--------------
 ;/localizations
 ;--------------
  (&localize f y ii tol-x tol-y maxit a b c sa sb sc g x u v  c w basta j plan)
  (&bind-too c u v w sa sb sc plan)
 ;---------------------------
 ;/make `f' hold the function
 ;---------------------------
  (if (symbolp f)		      
      (setf (symbol-function 'f) (symbol-function f))	; then clause
      (setf (symbol-function 'f) f))			; else clause
 ;----------------------------------------
 ;/get bounds for initial interval
 ;----------------------------------------
  (varbind 
   a (car ii)
   b (cadr ii))
  (unless (< a b) 
    (error "wrongly specified initial interval to `solve-monotone'"))
 ;--------------------------------------------------
 ;/define g(x) := f(x) - y , we will seek a zero of g 
 ;--------------------------------------------------
  (fbind g(x) (- (f x) y))
 ;---------------
 ;/define `basta'
 ;---------------
  (fbind basta ()
    (setq c (/ (+ a b) 2.0))
    (if (> (- b a) tol-x)
      (error 
       "solve-monotone: %s might be a solution if tol-x were %0.1g"
       c (- b a)))
    (if (> (o abs g c) tol-y)
      (error 
       "solve-monotone: %s might be a solution if tol-y were %0.1g"
       c (o abs g a)))
    (return-from solve-monotone c))
 ;-----------------------
 ;/begin iterative search
 ;-----------------------
  (psetq u (g a) v (g b))
  (loop 
   for j from 0
   if (> j maxit) do 
      (error "more than %d tries in `solve-monotone'" maxit)

   do (psetq sa (sgn u) sb (sgn v))
  ;----------------------------------------------
  ;/return if either a or b is already a solution
  ;----------------------------------------------
   if (= 0 sa) return a
   if (= 0 sb) return b
  ;-----------------------------
  ;/call `basta' if b-a <= tol-x
  ;-----------------------------
   if (<= (- b a) tol-x) do (basta)
  ;-----------------------------------------------------
  ;/determine whether root is within (a b) or outside it
  ;-----------------------------------------------------
   if (= sa sb) do (setq plan 'caminar) else do (setq plan 'narrow) end
  ;-------------------------
  ;/execute plan accordingly
  ;-------------------------
   do
   (case plan
    ;---------------------------------
    ;/case where root is outside (a b)
    ;---------------------------------
     ((caminar)
      (cond
	((< u v) (setq plan (if (> sa 0) 'left  'right)))
	((> u v) (setq plan (if (< sa 0) 'left  'right)))
	(t (basta)))
      (case plan
	((left)
  	   (psetq 
	    a (- (* 3 a) (* 2 b))
	    b a)
	   (psetq 
	    u (g a)
	    v u))
	((right)
	   (psetq 
	    a b
	    b (- (* 3 b) (* 2 a)))
	   (psetq 
	    u v
	    v (g b)))))
    ;--------------------------------
    ;/case where root is within (a b)
    ;--------------------------------
     ((narrow)
      (setq c (/ (+ a b) 2.0))
      (if (or (= a c) (= c b)) (basta))
      (setq 
       w (g c)
       sc (sgn w))
      (if (= 0 sc) (return c))
      (cond
       ((/= sa sc) (setq b c v (g b)))
       (t          (setq a c u (g a))))))))
 ;
 ; Possible improvements
 ;
 ;  add as keyword argts :lbound and :ubound between which to CONFINE search
 ;
 ;  so far this solves y=f(x) for x where f is a MONOTONIC
 ;  function.  It can easily be adapted to find a root of ANY function f given
 ;  an interval in which f changes sign. In fact implementing :lbound and
 ;  :ubound would take care of this.
 ;
 ; Notes
 ;
 ; The function `basta' takes over when iteration stops without reaching a
 ; solution.  It  issues an error mesage. 
 ;
 ; A problem that was cured was that for b-a very small c can equal one of
 ; them, then you just keep repeating the same interval (a b), now we catch 
 ; that and call basta.
 ;
 ; We will have trouble if f has a big flat stretch, but there is probably no
 ; point in trying to deal with this for now (we could do it by giving extra
 ; information to know whether it is monotone increasing or decreasing, or a
 ; memory of this) 
 ;
 ; If specify an initial b-a smaller than tol-x, it will stop immediately.
 ;
 ; c = midpoint of interval [a b] , w = g(c)
 ;
 ; The search strategy: 
 ;   if root is within interval, buscalo recursvely in the half where it is;
 ;   if it is outside, look in the adjacent interval of twice the size.
 ;
 ; Localization of variables can be important, in earlier version,
 ; `x' was left out, causing trouble. 
 ;
 ; The code in the first few lines is to allow the function-argument to be a
 ; symbol for the function as well as the fcn itself.  Another way would be to
 ; invoke the fcn using `funcall'.  Would this work/be better?  It probably
 ; would make the compiler happier, at any rate.


;::----------------------------------------------- 

;:: Functions devoted to printing/formatting 

(defun list-plist (L)
  "\
 Nicely prints the entries in the plist L, one pair per line. \
  "
  (cond
   (L 
    (tcl-or-elisp
     (princ (format nil "~s ~s ~%" (car L) (cadr L))) 
     (princ (format     "%s %s \n" (car L) (cadr L))))
    (list-plist (cddr L)))
   (t  "   ")))
 
(defalias 'print-plist 'list-plist)


(when *elisp*
(defun print-to-file (x file &optional overwrite)
  "\
 Prints any object to a file as if by `princ'.  If overwrite is nonnil
 it will overwrite the file, rather than appending to it.  The object is
 printed on its own line.
  "
  (write-region (format "\n%s\n" x) nil file (not overwrite)))   )
;
; In principle we want this for TCL too, someday


;; Some of the following aliases have misleading names and should be retired!
;; Also, some of the following functions are themselves largely useless (and
;; would be even more useless if elisp and TCL used the same `format'
;; conventions) 
;; Note: some builtins for printing are: princ prin1 print

(defun print-as-line (arg) 
  "\
 Prints a single lisp object ``for humans to read'', and preceded -- but not
 followed -- by a newline.  \
  "  
  (terpri) (princ arg))
 
;;; (defalias 'print-line 'print-as-line)
;;; (defalias 'printline  'print-as-line)

(defun print-list-as-line (&rest args)
  "\
 Prints multiple arguments as a single line ``for humans to read''.  
 Each argument is preceded by a space and the whole lot is followed by a 
 single newline.  Returns nil.
  " 
  (loop for $ in args do (princ " ")(princ $) finally (terpri)))
 
;;; (defalias 'print-list 'print-list-as-line) ; misleading name


(defun print-list-in-lines (&rest args) 
  "\
 Prints multiple arguments alternating with newlines 
 Produces output ``for humans to read'' (i.e. using `princ').
 The pattern is:  
                   <n><arg><n><arg>...<n><arg><n> 
 Returns nil. 
  "
  (loop 
    for $ in args 
    initially (terpri)
    do 
    (princ $)
    (terpri)))
 
(defalias 'print-lines 'print-list-in-lines)

(defalias 'breakline 'terpri)

