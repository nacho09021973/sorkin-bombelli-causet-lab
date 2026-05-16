;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;                                                          ;;;;;;;;
;;;;;;  All files in this directory are copyright 1997, 1998,   ;;;;;;;;
;;;;;;  1999 by Rafael D. Sorkin.  All rights reserved.         ;;;;;;;;
;;;;;;                                                          ;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;                        Time-stamp:<Dec 08 1998 14:14:07 (13933 31359)>
;
; The things in this file are, in one way or another, kept here for the sake of
; compatibility with TCL.  For example they may be needed by elisp but cause
; grief to TCL (in one of its incarnations), or they may be aliases needed to
; reconcile differences in function names between elisp and TCL, or it may be
; that they just are not needed by TCL and can conveniently be kept here.
;
; Several of these were once in `preparations.el' but were moved here to avoid
; the risk of the fcns using them never getting compiled properly.
;
; We probably want to load this file before bibliotek.macros, so don't use
; herein anything from there!  Or should we load it afterward so that we can in
; fact use macros from there?

;: Index
;|======================================================================
;| Functions, aliases and macros defined
;|
;|  A whole slew of simple aliases: 
;|     Assoc, Get, Member, Delete, Sort, Defun, 
;|     Compile, Let, find-symbol, defconstant
;| 
;|  Round        (needed since elisp lacks multiple return values)
;|
;|  Random
;|  random       (only temporary)
;|  random-from-elisp
;|  el:random    (to distinguish from one in Cl-package)
;|
;|  Read-from-string (needed since elisp lacks multiple return values)
;|
;|  (Once we also had things here to compensate for the lack of a proper    )
;|  (bacquote in elisp.  They are obsolete since emacs 19.29 and have gone  ) 
;|  (to the directory { old.code/backquote.stuff }.  They were called:      )
;|  (   Defmacro  Defsetf  el:bqf  el:bqm  el:bqh                           )
;|==========================================================================

;: Check we are really in elisp

 (unless *elisp* (error " Use { bibliotek.elisp } only with elisp"))


;: A bunch of simple aliases

; The following aliases are meant to allow elisp code to be compatible with
; TCL (true common lisp) conventions.

(defalias 'Assoc  'assoc*)
(defalias 'Get    'get*)
(defalias 'Member 'member*)
(defalias 'Delete 'delete*)  
(defalias 'Sort   'sort*)
(defalias 'Defun  'defun*)
(defalias 'Compile 'byte-compile)
(defalias 'Let    'lexical-let)
(defalias 'find-symbol 'intern-soft)

(defalias 'defconstant 'defconst)	
 ; For some strange reason, this alias has problems when being loaded after
 ; compilation.

;: Some further aliases, macros and functions

(defun Round (x y)
  "\
 Like the TCL version, except returns only first return value, in other
 words, (Round x y) => the nearest integer to x / y"
  (car (round* x y)))

;:: Rearranging the names of random number generators for consistency with TCL

; We use `Random' for the cl-package version that acts like TCL's `random'.
; The version built into elisp, we rename to `el:random'.
; We also define `random-from-elisp' to mimic TCL's random using the elisp
; built-in.
; The name `random' we don't use at all, since using it for the CL version has
; the problem that if this file is re-executed, then el:random will become the
; wrong thing; this way, we'll just get an error.

(defalias 'Random 'random*)

(defalias 'el:random  (symbol-function 'random))

(defun random (&optional MMM)
  " This is now a dummy function. Use `Random' or `random-from-elisp' instead."
  (error "Don't use `random', use `Random' or `random-from-elisp' instead."))

(defun random-from-elisp (limit)
  "\
 A random number in [0 limit).
 If `limit' is integral then so is the value returned.
 If `limit' is floating point then the possible values are
 are actually discrete, with a spacing of one trillionth of
 the size of the interval.
 To initialize the seed do (el:random t).
 Beware: gives garbage if `limit' is negative or zero.
 This is concocted from elisp's `random' which returns integers only 
 (they can be negative), and which gives garbage if its argument is float!!
  "
  (unless (= 1000000000000 1000000000000.0)
    (error "The integer 1000000000000 is too big for this machine!")) 
  (cond
   ((integerp limit) (el:random limit))
   (t
    (* (/ limit 1000000000000.0) (el:random 1000000000000)))))


(defun Read-from-string (string) 
  " Invokes the corresponding elisp fcn and discards the second return value."
  (car (read-from-string string)))

