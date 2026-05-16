;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;                                                          ;;;;;;;;
;;;;;;  All files in this directory are copyright 1997, 1998    ;;;;;;;;
;;;;;;  by Rafael D. Sorkin.     All rights reserved.           ;;;;;;;;
;;;;;;                                                          ;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;:==============================================================
;: Implement patches to elisp, cl, etcetera if they are needed,
;: i.e. if not yet done systemically
;:==============================================================
                             ; Time-stamp:<Sun Dec 22 17:23:46 (12989 49922)>

;--------------------------------------------------------------------
;: Bug affecting lexical-let and (thereby) labels as of emacs 19.34.1
;--------------------------------------------------------------------

;::  Define Labels to be flet

; This was originally to evade a bug in cl-19
; which in particular affects my macro `deff'
; Howver, it might be a good idea in general to replace `labels' with `flet'
; in elisp, since we don't rely on lexical scoping anyhow, and `flet' is
; probably faster than `labels'
; To make `Labels' really be `labels' just use the alternative defalias
; provided: 
 
(defalias 'Labels 'flet) 
;; (defalias 'Labels 'labels) 

(message " aliased `Labels' to `flet'")
 
;::  And now we provide Dave Gillespie's actual patch for lexical-let

(defmacro lexical-let (bindings &rest body)
  "(lexical-let BINDINGS BODY...): like `let', but lexically scoped.
The main visible difference is that lambdas inside BODY will create
lexical closures as in Common Lisp."
  (let* ((cl-closure-vars cl-closure-vars)
	 (vars (mapcar (function
			(lambda (x)
			  (or (consp x) (setq x (list x)))
			  (cl-push (gensym (format "--%s--" (car x)))
				   cl-closure-vars)
 			  (set (car cl-closure-vars) [bad-lexical-ref]) ;; NEW
			  (list (car x) (cadr x) (car cl-closure-vars))))
		       bindings))
	 (ebody 
	  (cl-macroexpand-all
	   (cons 'progn body)
	   (nconc (mapcar (function (lambda (x)
				      (list (symbol-name (car x))
					    (list 'symbol-value (caddr x))
					    t))) vars)
		  (list '(defun . cl-defun-expander))
		  cl-macro-environment))))
    (if (not (get (car (last cl-closure-vars)) 'used))
	(list 'let (mapcar (function (lambda (x)
				       (list (caddr x) (cadr x)))) vars)
	      (sublis (mapcar (function (lambda (x)
					  (cons (caddr x)
						(list 'quote (caddr x)))))
			      vars)
		      ebody))
      (list 'let (mapcar (function (lambda (x)
				     (list (caddr x)
					   (list 'make-symbol
						 (format "--%s--" (car x))))))
			 vars)
	    (apply 'append '(setf)
		   (mapcar (function
			    (lambda (x)
			      (list (list 'symbol-value (caddr x)) (cadr x))))
			   vars))
	    ebody))))

(message " patch applied to `lexical-let'")

;-------------------------------------------------------------
;: Bug in equalp  as of 19.30 fixed by 19.34
;-------------------------------------------------------------
;   For the bug in `equalp' (CL package) present as of emacs 19.30 
;   From: Richard Stallman <rms@gnu.ai.mit.edu>
;   Subject: Re: bug in CL package
;   Thanks.  Please try this fix.
;-------------------------------------------------------------

; apply patch only if needed

(if 
    (equalp (list 0) (list 1))

(progn
(defun equalp (x y)
  "T if two Lisp objects have similar structures and contents.
This is like `equal', except that it accepts numerically equal
numbers of different types (float vs. integer), and also compares
strings case-insensitively."
  (cond ((eq x y) t)
	((stringp x)
	 (and (stringp y) (= (length x) (length y))
	      (or (equal x y)
		  (equal (downcase x) (downcase y)))))   ; lazy but simple!
	((numberp x)
	 (and (numberp y) (= x y)))
	((consp x)
;;	 (while (and (consp x) (consp y) (equalp (cl-pop x) (cl-pop y))))
   	 (while (and (consp x) (consp y) (equalp (car x) (car y)))
   	   (setq x (cdr x) y (cdr y)))
	 (and (not (consp x)) (equalp x y)))
	((vectorp x)
	 (and (vectorp y) (= (length x) (length y))
	      (let ((i (length x)))
		(while (and (>= (setq i (1- i)) 0)
			    (equalp (aref x i) (aref y i))))
		(< i 0))))
	(t (equal x y))))  

(message " Patch applied to `equalp' ")))
