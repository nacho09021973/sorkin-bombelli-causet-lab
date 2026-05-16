;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;                                                          ;;;;;;;;
;;;;;;  All files in this directory are copyright 1997, 1998,   ;;;;;;;;
;;;;;;  1999 by Rafael D. Sorkin.  All rights reserved.         ;;;;;;;;
;;;;;;                                                          ;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;==========================================================================
;                  General Preparations for elisp 
;==========================================================================
;
; (Hay tambien counterparts like `preparations.gcl' for use with dialects of
;  TrueCommonLisp) 
;                      ; Time-stamp:<Dec 21 1998 01:42:56 (13949 60912)> 
;=======================================================================

;: Decide whether or to recompile everything 

(defvar recompile nil)
; (setq recompile  t )
; (setq recompile nil)

(message "================================\nBeginning preparations for elisp")

;=======================================================================

;: Table of Contents 
 ;  
 ;  Define some global constants
 ;  Declare and initialize some global variables     
 ;  Adjust various limits   
 ;  Kwa ajili ya TCL
 ;  Try to define Infinity% and NaN%
 ;  Adjust Load Path
 ;  Load various libraries (with recompilation if requested)
 ;  Other 
;----------------------------------
;: Define some global constants 
;----------------------------------
(progn
  (defconst *lisp-type* 'elisp)		; elisp with CL extensions
  (defconst *elisp* t)
  (defconst *TCL* nil)			; TCL = True Common Lisp
  (defconst *GCL* nil)			; GCL = Gnu  Common Lisp
  (defconst *gcl* nil)			; gcl = Gnu  Common Lisp
  (defconst *clisp* nil)
  (defconst *sun-lisp* nil)		; sun version of lucid lisp
  "some constants defined")
;----------------------------------------------
;: Declare and initialize some global variables 
;----------------------------------------------
(defvar *carefully* t 
  "\
 When true, makes various functions operate with more checking:
 poset functions may check for transitivity or sorting, matrix 
 multiplication may check that dimensions of factors match, etcetera." )
 ; 
 ; Note: the default value (t) takes effect only if *carefully* is void.

(defvar *mrv* nil
  " Place to deposit multiple return values")

;-----------------------
;: Adjust various limits 
;-----------------------
(setq 
 max-lisp-eval-depth (expt 2 14)   ; (was 2^14) recursion depth 
 max-specpdl-size    (expt 2 16)   ; (was 2^16) number of variable bindings 
 gc-cons-threshold   (expt 2 24)   ; frequency of garbage collection (in bytes)
 print-length nil		   ; truncates printing if nonnil
 print-level  nil)		   ; limits nesting depth in printing

(message "  various limits increased governing recursion depth, etc.")  

 ; The reduced limits of 2^13 were for umoja-suhep
 ; They were best possible as of Thu Dec 25 1997
 ; but doing ulimit -s 32000 seemed to cure the problem
 ;
 ; Reducing the frequency of garbage collection, gains perhaps 
 ; a factor of 2 to 3 in speed. 
 ;
 ; Newer emacs versions admit integers up to 2^27-1, at
 ; least.  (Older ones could be limited to 2^23-1.)

;--------------------------------------------------------------
;: Kwa ajili ya TCL 
;--------------------------------------------------------------

(defun in-package (&optional $) 
  "A dummy function as far as elisp is concerned" nil)

 ; More things related to TCL compatibility are in `bibliotek.elisp.el'

;-------------------------------------------------
;: Try to define Infinity% and NaN% 
;-------------------------------------------------
(condition-case 
 error-message
 (progn
  (require 'cl)
  (defconst Infinity% 1e5000    "The ``infinite number''")
  (defconst Infinity  Infinity% "The ``infinite number''")
  (defconst NaN% (- Infinity% Infinity%) "``Not a Number''")
  (defconst NaN NaN%)
; [Now test for them really being +Infinity and NaN]
  (unless 
    (and
     (numberp NaN)
     (numberp Infinity)
     (not(or
       (> NaN 0)
       (< NaN 0)
       (= NaN 0)))
       (= Infinity (+ Infinity Infinity))
       (/= Infinity 0)
       (> Infinity 0))
    (message
     "  Infinity and NaN came out WRONG:\n   Infinity => %s and NaN => %s." 
     Infinity% NaN%)
    (message "   Both will be made symbols (rather than IEEE floats)")
    (setq 
     Infinity% 'Inf
     NaN%      'NaN
     Infinity   Infinity%
     NaN        NaN%) 
    "  Have not defined Infinity% or Nan% properly"))
 (error
  (message 
   "  failed to define Infinity%% or NaN%%, error message was:\n  %s." 
   error-message)))
  ;  
  ; Notes
  ;
  ; Except for Infinity% and NaN%, all other numerical constants are loaded
  ; from the file "bibliotek.constants.el" 
  ;
  ; Apparently, one NaN need not be equal (or = or eql or equalp) to another!
  ;
  ; We trap errors here because the attempt to define Infinity etc sometimes
  ; produces arithmetic errors.  The checks are done because on some machines
  ; you get actual numbers rather than the IEEE things.  
  ; One implementation that doesn't do what we want is DEC emacs.
 
;-------------------------------------------------------------------------
;: Adjust Load Path 
;-------------------------------------------------------------------------

(if (not(boundp '*nyumbani*))
 (progn    
   (setq 
     *nyumbani*
     (cond
      ((equal system-name  "umoja.syr.edu") "/home1/sorkin/")
      ( t                                   "~/")))
   (setq
     load-path 
     (append
      load-path
      (list (concat *nyumbani* "lisp")))))) 

 ; This makes sure that our lisp directory is in the load path, and any other
 ; directories we might need can also be added (e.g. "lisp/cl-source.files").
 ; Actually all this is done already by the .emacs file, however if we invoke
 ; emacs using `emacs -q' then `.emacs' doesn't get loaded, and the above
 ; kicks in then.  
 ; The variable `*nyumbani*' should be where we want our
 ; effective base directory to be.
 ;
 ; In the load-path `nil' represents the current directory.  We do NOT include
 ; it because it leads to trouble when there is a local file sharing a name
 ; with one the system is looking for (e.g. "outline" which can be confused
 ; with the system's "outline.elc").  Of course this entails the opposite
 ; danger of a system file pre-empting a local one.  To avoid that, either use
 ; more complete pathnames for local files, or give them names the system
 ; wouldn't be likely to employ.
 ;
 ; There are many ways to find out which machine you are on without executing a
 ; a shell command, or using `getenv', eg the variable `system-name' has it.  
 ; To obtain the version number of emacs use, eg
 ;      (read (substring emacs-version 0 5))

;------------------------------------------------------------
;: Load various libraries (with recompilation if requested) 
;------------------------------------------------------------
; Notes
;
;   Obviously the order of loading can be important.
;
;   You can use M-x byte-force-recompile to recompile all the .elc files
;   (This exists as of emacs version 19.30)
;   But it acts on ALL SUBDIRECTORIES, which we might not want, also it seems
;   betters to control the order of compilation.
;
;   Most lisps use plain `load', and don't recognize `load-file'.
;
;--------------------------------------------
;:: Load the CL package that comes with emacs
;--------------------------------------------
(progn
(load "cl")				; master file (calls others if needed)
(load "cl-macs")			; macros
(load "cl-seq")				; concerning sequences 
(load "cl-extra")			; larger, more complex functions
(require 'cl)
(require 'cl-19))

 ; It is better to make sure all the CL files are loaded at this stage. 
 ; See the warning about redefining functions from cl-package. 
 ; Must make sure all are byte-compiled, if want speed.

;------------------------------------------
;:: Load a file of patches for elisp and cl
;------------------------------------------

(load "bibliotek.elisp.patch")

  ; Note: this will NOT be compiled below, compile it by hand when it is
  ; changed!

;-------------------------------------------------------
;:: Load my biblioteks (with recompilation if requested)
;-------------------------------------------------------
;
(when recompile

  (load  "bibliotek.elisp.el")
  (load  "bibliotek.macros.el")
  (load  "bibliotek.constants.el")	; We do NOT compile this one
  (load  "bibliotek.general.el")
  (load  "bibliotek.poset.el")
  (load  "bibliotek.float.el")
  (load  "bibliotek.extras.el")

  (byte-compile-file "bibliotek.elisp.el"   'load)
  (byte-compile-file "bibliotek.macros.el"  'load)
  (byte-compile-file "bibliotek.general.el" 'load)
  (byte-compile-file "bibliotek.poset.el"   'load)
  (byte-compile-file "bibliotek.float.el"   'load)
  (byte-compile-file "bibliotek.extras.el"  'load)

  (message "Bibliotek's recompiled and loaded"))

 ; The nonnil second argument in the above means to load after compiling

(unless recompile

  (load  "bibliotek.elisp")  
  (load  "bibliotek.macros")
  (load  "bibliotek.constants.el")	; We do NOT compile this one
  (load  "bibliotek.general")
  (load  "bibliotek.poset")
  (load  "bibliotek.float")
  (load  "bibliotek.extras")

  (message "Bibliotek's loaded without recompilation"))

(setq recompile nil)

;--------------------
;: Other 
;--------------------

;; (setq  *gensym-counter* 0)  
;;
;;  Uncomment this if you want to make cl-19 behavior repeatable for debugging
;;  purposes. (It's not good to do otherwise, since it carries a risk of name
;;  collisions.)

(provide 'preparations)	    ; Here "provide" really means "I have provided"

(message "Preparations for elisp completed\n================================")
 
;:                   - end -
