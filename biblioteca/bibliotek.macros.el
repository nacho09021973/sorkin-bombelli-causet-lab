;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;                                                          ;;;;;;;;
;;;;;;  All files in this directory are copyright 1997, 1998    ;;;;;;;;
;;;;;;  by Rafael D. Sorkin.     All rights reserved.           ;;;;;;;;
;;;;;;                                                          ;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; Time-stamp:<17:22:19 Jul 11 1998 (13735 55179)>

;: NEEDED CHANGES ARE FLAGGED WITH `;;;'

;: This is the file { bibliotek.macros.el }
;
; Macros are segregated here because they are needed by many other
; biblioteks, hence are best loaded first.
;
; For tests of the macros defined herein see "tests.and.timings/bib.macros.el"

;: CAUTION Don't invoke herein fcns that are defined in other biblioteks!

;: Advice on compiling biblioteks
;
;  In compiling biblioteks from scratch, by far the best policy seems to be to
;  LOAD the files FIRST and only THEN COMPILE them.
;  If you compile without first loading, the functions and macros seem not to
;  know about each other, even about ones occurring earlier in the file. 
;  Probably this can be overcome by use of `eval-when' (there would be no harm
;  in wrapping everything in it, in fact), but experience shows that it's far
;  easier and better just to load before compiling.

;: Advertencia:  Compiling `deff' and `kwa' is SLOW in clisp and elisp

;: Roster of macros and functions defined (not in order)
;| 
;|  setqq       (like setq but second args are also not evalled)
;| 
;|  compose     (composition of functions acting on SINGLE arg)
;|
;|  image      (a more convenient syntax for using `map', `mapcar')
;|
;|  kwa        (a simple looping macro effectively superseded by `loop')
;|
;|  arith-to-lisp/atl   (Allows the use of normal infix syntax in arithmetic
;|                       expressions) 
;|
;|  The csynf family 
;|
;|    csynf + helpers (csynf = clearer-syntax-facilitator)
;|
;|    deff            (form of defun using csynf syntax)
;|
;|    toggle-csynf-debugging (helps in debugging functions defined with `deff')
;|
;|  localize-using-name/LUN   (macro used to make variables
;|                            local by creating new symbols for them.)
;|  replace-using-name        (used in the above)
;|
;|  tcl-or-elisp (to use different code depending on TCL or elisp)
;|
;|  defconst-f    (like defconst, but a function not a macro)
;|                (being developed elsewhere?)                  
;|
;|  Time and time
;|
;|    Time   (does garbage collection first, then `time')
;|    time    (a very crude timer, giving only elapsed time)
;|    get-time (get current time in seconds)
;|    pp-elapsed-time
;|
;|  retrieve-error-messages-from  (no longer needed)
;|
;|  Some aliases which are defined in 
;|     preparations.el or patches.elisp or bibliotek.elisp.el 
;|
;|    Defun    (an alias for defun* in CL extension)
;|    Let      (an alias for lexical-let in CL extension)
;|
;|    Labels   (for now it's an alias for flet, since labels has trouble)  
;|  
;|========================================================================

(in-package 'user)

;:----------------------------
;: The macro `setqq'  
;----------------------
(unless *clisp*

(defmacro setqq (&rest plist)
  " Like setq but *neither* element of each argument pair is evaluated"
  (cons 
   'progn
   (loop 
    for $ on plist by 'cddr
    collect `(setq ,(car $) (quote ,(cadr $))))))   )
   ;
   ; Note: we don't use `set' since in TCL it will act on the GLOBAL
   ; binding of symbol, whereas `setq' acts on the local binding

;; stripped down version for CLISP (since lacks loop)
(when *clisp*
(defmacro setqq (symbol object)
  " Like setq but *neither* arg is evaluated.  Takes only one arg pair"
  `(setq ,symbol (quote ,object))) )
 
;:----------------------------
;: The macro `tcl-or-elisp'
;------------------------- 
(defmacro tcl-or-elisp (tform eform)
  "Return the first form if in TCL, the second if in elisp"
  `(cond
    (*TCL*   ,tform)
    (*elisp* ,eform)))
;:----------------------------
;: The macro `compose'
;-------------------------
(defmacro compose (&rest H)
  " 
 A macro to write more simply the composition of a string of functions 
 acting on a single argument, eg: `(compose sin cos tan sqrt 5)' or,
 using the convenient alias `o':        `(o sin cos tan sqrt 5)'.
  "       
  (if (> 2 (length H)) (error "Error: too few arguments to `compose'."))
  (Labels 
   ((inner (H)
      (cond 
       ((= 2 (length H)) H)
       (t (list 
           (car H)
           (inner(cdr H))))))) 
   (inner H)))

(defalias 'o 'compose)    
 
;:----------------------------------------------------------------------
;: Macros for localizing symbols used within functions (LUN and friends)
;-----------------------------------------------------------------------
;:: Notes on the following macros and functions
;
; A better method of localization MAY be developed someday in
; ~/lisp/developing/new.localization/  See that file for ideas.
;
; `replace-using-name' uses a name given to it to deterministically construct
; new symbol names from the old ones.  This is in one sense less safe than
; using `gentemp' since it doesn't check for each name being truly new,
; but its big advantage is that it should work better with compiled functions,
; since when they are reloaded, `gentemp' will NOT be re-executed.   
;
; `replace-using-name' doesn't really need to be separate from
; `localize-using-name' (especially if you use `sublis', which does
; simultaneous substitutions of many symbols)
;
; In `replace-using-name' we could use `make-symbol' rather than `intern' to
; get uninterned symbols.  These would be "even more local" -- before
; compilation and loading.  Is there any harm in it?  Well, maybe a false
; impression of security.
;
; The "temps" pair below probably is never used, and could be trashed.
;  For now just comment it out and see what happens

(defun replace-using-name (name sym-list vitu)
  "\
 Arguments are:  name (a string),  sym-list (a list of symbols),  vitu .
 We return 
 vitu unchanged except that occurrences of the elements of sym-list
 are replaced with other symbols whose names are formed from those of the
 originals using name.  (It is intended that the new name be unique.  If by
 accident it is not, then the existing symbol with that name becomes the
 replacement.)    \
  "
  (if (null sym-list) vitu
    (subst 
     (intern
      (concatenate 
       'string    "-<"  (symbol-name (car sym-list))  ">-@-"  name)) 
     (car sym-list) 
     (replace-using-name name (cdr sym-list) vitu))))  
  ; 
  ; Other forms which look perhaps nicer (current is  -<x>-@-t-close ):
  ;
  ; -<x>-%-{t-close}
  ; -<x>-@-{t-close}
  ;
  ;  We do NOT localize the variables in `replace-using-name', since it appears
  ;  to be unnecessary.  If wrong then just do it by hand, either changing
  ;  eg "vitu" to "-<vitu>-@-replace-using-name" etc., or wrapping the whole
  ;  thing in a `sublis'.)  
  ;  Symbols to be localized would be: name sym-list vitu


(defmacro localize-using-name (-<name>-@-{LUN} 
                               -<syms>-@-{LUN} 
                        &rest -<forms>-@-{LUN})
  "\
 Replaces the symbols with new ones whose names are made from the
 given name (purpose being to make them effectively local).
 Example: 
            (localize-using-name \"t-close\" (x y j k) FORMS) 
 CAUTION: 
 This replaces ALL instances of the symbol.  So be sure not to give your 
 local symbols names which collide with ``external'' (or ``special'' or
 whatever) variables which you want to refer to in the function body, or with
 symbols you want to use in a quoted context (as plist indicators for
 example); finally of course, don't localize symbols which are names of
 external functions you want to use, like `sqrt'. 
  "
  (cons 
   'progn 
   (replace-using-name -<name>-@-{LUN} -<syms>-@-{LUN} -<forms>-@-{LUN})))
 
(defalias 'LUN 'localize-using-name)
 
;:------------------------------------------------------
;: The csynf family 
;;---------------
 ;
;:: SOME COMMENTS on this family
 ;
 ; Further testing may be needed!
 ;
 ; The name csynf stands for "clearer-syntax-facilitator"
 ;
 ; Method is described on sheet of paper somewhere, and below incompletely
 ;
 ; For TCL should we use flet or labels?   The diff seems to be that labels
 ; allows the new fcn to be recursive, whereas flet wouldn't.  So for now
 ; we choose to use `labels' !  In either case the newly defined fcns cannot
 ; refer to each other mutually, which is perhaps nuisance.  (They could if we
 ; allowed multiple pairs within fbind (like with varbind) and used `labels'
 ; as the binder)
 ;
 ; csynf-inner is kept separate form csynf for clarity.
 ; auxiliary fcns also kept separate for now, since seems to speed things up
 ; (at least it does so a lot before compilation) (why??)
 ; if put them inside, be sure to use `labels' to allow recursion
 ;
 ; Neither `varbind' nor `fbind' are more than "pseudo-functions" of course.
 ; Below, we define them to give errors if you really try to evaluate a
 ; pseudoform containing them.  We make them be macros rather than functions,
 ; so that the error messages will be clear.  (If they were functions, errors
 ; could arise earlier from trying to eval their args.) 
 ;
 ;
;;: POSSIBLE IMPROVEMENTS
 ;
 ;  Make `&bind-too' work even if it is not at top level, so it can be used
 ;  within an `fbind' for example.
 ;
 ;  Should we add a `varbind-parallel' which uses `let' rather than `let*'?
 ;
 ;  In addition to (&localize...) we could have:
 ;   (&doc ....) to allow more flexibility in including doc-string
 ;   (&aux ....) but this would seem to be redundant now, since we have
 ;   `varbind' 
 ;
 ;  Make also a `deffmacro'  to allow same stuff within macro defs?
 ;  But probably easier just to wrap a `csynf' around them !
 ;
 ;
;;: EXPLANATIONS (The functions/macros and what they do)
 ;
 ; `csynf-binder-p' tests whether its arg is a list beginning 
 ;                  with `varbind' or `fbind' 
 ;
 ; `csynf-process-vbind-list' puts parens around the pairs 
 ;                  for feeding to `let*'
 ;
 ; `csynf-inner' is the guts of the macro `csynf'
 ;
 ; Summary of how it all works
 ;
 ;  Let X Y Z denote any lisp objects, with X* Y* Z* being successions
 ;            (possibly empty) of such.
 ;
 ;  Let A B C denote any lisp objects EXCEPT for lists beginning with one of
 ;            the "binders" `varbind' or `fbind', and let  A* B* C* again
 ;            denote successions of such.
 ;
 ;  For example "(setq A*)" stands for (setq A1 A2 A3...)
 ;
 ;  We want something like 
 ;
 ;   (A* (fbind B*) C*)   to become   (A* (labels ((B*)) C*))  
 ;
 ; or rather to become this with A* etc also processed.  To implement this we
 ; take the list H = (X Y Z ...) and examine X.  If X is anything but a list
 ; beginning with `fbind'or `varbind' we just pass over it, or rather we
 ; process it recursively and put it back in place.  (Exception: if it is a
 ; list beginning with `quote' we don't process it.)
 ; When we come to a list which does begin with a "binder" then we treat it as
 ; illustrated in one of these special cases:
 ;
 ;  ((fbind X*) Y*)  -->  ((labels ((X*)) Y*)) 
 ;
 ;  ((varbind X Y ...) Z*)  --> ((let* ((X Y)...) Z*))
 ;  
 ; (In these examples the recursive processing of X Y etc is not indicated
 ; explicitly).
 ;
;;: Notes on `deff'
 ;
 ; `deff' is a macro for incorporating "clearer syntax" and
 ; `localize-using-name' into the function definition itself, so they don't
 ; have to be invoked explicitly 
 ;
 ; In defining it, we make the `csynf' come *before* the `defmacro' so that the
 ; effect of any `(interactive)' (in elisp) won't get messed up.
 ;
 ; The vble `loc-form' is the list (&localize sym1 sym2 ...)
 ; process-loc-list should yield (sym1 nil sym2 nil...)
 ; loc-list will be nil when (&localize sym1 sym2 ...) is absent
 ;
 ; we could also allow a form (&doc...) for the documentation, but
 ; what syntax to use within it?  if have it at all might as well allow
 ; multiple strings and also ``...'' to become \"...\"
 

(defun csynf-binder-p (A) (and (listp A) (memq (car A) '(varbind fbind))))

(defun csynf-quote-p  (A) (and (listp A) (eq (car A) 'quote)))

(Defun csynf-process-vbind-list (B) 
  (if 
    (null B) 
    (return-from csynf-process-vbind-list nil))
  (setq B (csynf-inner B))           ; to allow `varbind' etc within varbind
  (cons
   (subseq B 0 2)
   (csynf-process-vbind-list (cddr B))))

(Defun csynf-inner (H)
  " This function (not macro) is the guts of the macro `csynf'.  It is kept
 separate from the latter only for legibility.  (This separation would not be
 necessary if we could use the clearer `fbind'-syntax within csynf itself!)" 
 ;-------------------------------
 ; trivial, recursion-ending case
 ;-------------------------------
  (when (null H) (return-from csynf-inner nil))
 ;------------------------------------------------------------
 ; case where first form in H is `(quote...)' or
 ; not a list at all : don't touch it
 ;------------------------------------------------------------
  (when 
      (or
       (o not listp car H)
       (o csynf-quote-p car H))
    (return-from csynf-inner    
      (cons 
       (car H) 
       (o csynf-inner cdr H))))
 ;-----------------------------------------------------------
 ; case where the first form in H is a list, but not an fbind 
 ;     varbind (or "quotation")
 ;-----------------------------------------------------------
  (when (o not csynf-binder-p car H)
    (return-from csynf-inner
      (cons 
       (o csynf-inner car H)
       (o csynf-inner cdr H))))
 ;-----------------------------------------------------------------------------
 ; remaining case: first form in H is list beginning with `fbind' or `varbind'
 ;-----------------------------------------------------------------------------
  (case (caar H)
   ;----------------
   ; subcase: fbind
   ;---------------
    ((fbind)
     `((Labels                          ; notice we don't use flet
         ( ,(o csynf-inner cdar H))
          ,@(o csynf-inner  cdr H))))
   ;------------------
   ; subcase: varbind 
   ;------------------
    ((varbind)                         
     `((let* 
          ,(o csynf-inner csynf-process-vbind-list cdar H)
         ,@(o csynf-inner                           cdr H))))))    

(defmacro csynf (&rest H)
  " A macro whose name means {c}learer{syn}tax{f}acilitator.  It allows you
 to use `varbind' and `fbind' instead of `let*' and `labels'.  They have 
 better names, demand fewer parentheses, and don't need to be ``wrapped
 around'' the forms in their scope.  Rather their scope extends until the 
 end of the form in which they occur.  (In this sense each form acts as 
 a ``block'', as far as csynf is concerned.)
 Usage: 
           (csynf FORM FORM...FORM)     (Value will be that of final form.)
 Examples:
             (varbind x 1 y 2 z x)   
             (fbind phi (x y) (+ x y))
 Notice here
 that `varbind' can bind multiple args, and does so sequentially, rather 
 than in parallel, whereas `fbind' binds only a single function.
 Beware: csynf can run into trouble when combined with backquote (so in
         particular in defining macros, see file :WARNINGS for advice)
 Comment: the reason we make `fbind'  use `labels' rather than `flet' is that 
 in TCL, functions defined by `flet' apparently can't be recursive.
  "
  (cons 'progn (csynf-inner H)))
 ;
 ; See further explanations above

(defalias 'clearer-syntax-block 'csynf)

(defmacro varbind (&rest args) 
  " This is just a dummy function so far.  If `clearer syntax' were built into
 TCL then this could be a real function replacing let*.  As it is, we just
 make it signal an error, since it should never be used as such (except
 for debugging purposes).  See `csynf' for more information."
  (declare (ignore args)) 
  (error "Error: `varbind' is not a real function or macro in elisp or TCL"))

(defmacro fbind (&rest args) 
  " This is just a dummy function so far.  If `clearer syntax' were built into
 TCL then this could be a real function replacing `flet' or rather `labels'.
 As it is, we just make it signal an error, since it should never be used as
 such (except for debugging purposes).  See `csynf' for more information."
  (declare (ignore args))
  (error "Error: `fbind' is not a real function in elisp or TCL"))
 
;:: Advertencia: Defining `deff' is  S L O O O W  in clisp and elisp

;------------------
;:: The macro `deff'
;------------------
(csynf
(defmacro deff (&rest H) 
  " This is an extension to `defun'.  You define a function as usual except 
 that
  (a) you can use the `varbind' and `fbind' syntax within it,
  (b) you can include a ``pseudoform'' looking like (&localize x y z...),
  (c) you can include -- at top level only! -- an (&bind-too ...) pseudoform. 
 If you include (b) the symbols x y z... will automatically be localized 
 by `localize-using-name', the name being that of the function itself.
 If you include (c) then that pseudoform (wherever it may be, as long as it
 is at top level) will turn into the equivalent of (varbind x nil y nil...).
 This last is useful mainly for preventing compiler warnings (and to that
 extent should not really be needed at all), but can be convenient on its own
 now and then.
    See the documentation of `localize-using-name' for more on what that macro
 does, including some cautions.  
    In the elisp case, `deff' also serves to incorporate the CL extensions 
 to elisp's `defun'.   
    The name ``deff'' stands for ``define-function'' or ``defun-fancy''."
 ;----------------------------------------
 ; Define auxiliary function for &bind-too
 ;----------------------------------------
  (fbind process-bind-list (B) 
    (cond
     ((null B) nil)
     (t (append
	 (list (car B) nil)
	 (process-bind-list (cdr B))))))
 ;-----------------------------------
 ; Record name of function being defined
 ;-----------------------------------
  (varbind name (o symbol-name car H))
 ;--------------------------------------------------
 ; Find the `&localize' pseudoform if there is one, extract the symbols to be
 ; localized and delete the pseudoform if it was there.
 ;--------------------------------------------------------------------------
  (varbind 
   loc-form (find '&localize H :key 'car-safe)
   loclist (cdr loc-form))
  (if loc-form (setq H (delq loc-form H)))
 ;--------------------------------------------------
 ; Find the `&bind-too' pseudoform if there is one
 ;--------------------------------------------------
  (varbind 
   bind-too-cons (Member '&bind-too H :key 'car-safe)
   bind-list nil)
 ;--------------------------------------------------------------------------
 ; extract the symbols to be varbound and convert the @ pseudoform to do so
 ;--------------------------------------------------------------------------
  (when bind-too-cons
    (setf 
     bind-list (cdar bind-too-cons)  
     (car bind-too-cons) (cons 'varbind (process-bind-list bind-list))))
 ;------------------------------------------------
 ; rebuild the function def with localized symbols
 ;------------------------------------------------
  `(localize-using-name ,name ,loclist
     (csynf 
      (Defun ,@H)))))
 ;
 ; See further explanations above

;----------------------------
;:: Nice feature for debugging
;----------------------------

(defvar *csynf-debug* nil "true when varbind redefined to setq etc")
(defvar *csynf-save-fbind* nil "holds original value of fbind")
(defvar *csynf-save-varbind* nil "holds original value of varbind")

(defun toggle-csynf-debugging () 
  " Temporarily bind `varbind' to `setq' and `fbind' to `deff'"
  (interactive)
  (cond 
   ((not *csynf-debug*)
     (setq *csynf-debug* t)
     (setf 
       *csynf-save-varbind*  
       (symbol-function 'varbind))
     (defalias 'varbind 'setq) 
     (setf 
       *csynf-save-fbind*  
       (symbol-function 'fbind))
     (defalias 'fbind 'deff)
     (message "csynf debugging turned ON"))
   (*csynf-debug*
    (setq *csynf-debug* nil)
    (setf 
      (symbol-function 'varbind)
      *csynf-save-varbind*)
    (setf 
      (symbol-function 'fbind)
      *csynf-save-fbind*)
    (message "csynf debugging turned OFF"))))
 
;:----------------------------------------------------------
;: the macro `image'  (checked out pretty well, not totally) 
;;----------------------------------------------------------
                 (LUN 
          "image" (args D D-supplied F F-supplied type type x)
(defmacro* image (&rest args
                 &aux 
                 D (D-supplied nil)
                 F (F-supplied nil)
                 (type 'list)
                 (x '$))
  "\
 Treats an expression F as a function of the symbolic argument `$' (or the
 arg provided) and applies it to each element of a Sequence D, returning the
 resulting sequence of values as a list (or whatever called for).  The clauses
 can be in any order (as for `loop') and many different ``keywords'' can be
 used to introduce them,
 for example
               (image of (* $ $) on '(1 2 3))
               (image :on L :of (cons $ $))
               (image for x in L of (* x x))
               (image of (+ $ 5) on M out-type vector)
               (image vector on L of (* $ $))
               (image of (+ $ 5) on M output nil)

 Here is a full list of the `keywords' accepted grouped into lists of 
 synonyms:
        (domain :domain dom :dom on :on in :in)
        (of :of expression :expression expr :expr form :form formula :formula )
        (out-type :out-type out-to :out-to output :output out :out)
        (argument :argument arg :arg parameter :parameter for :for)
 The possible 
 out-types are `vector' `list' and `sequence' (though TCL may not accept last).
 They can be specified on their own as well as in, eg, `type list'.\
 "
  (while args 
   (case (car args)
     ((on :on in :in dom :dom domain :domain)           ; add `from' `to'?
      (setq    
       D (cadr args)
       D-supplied t
       args (cddr args)))
     ((of :of form :form formula :formula expr :expr expression :expression)
      (setq    
       F (cadr args)
       F-supplied t
       args (cddr args)))
     ((out :out out-to :out-to out-type :out-type output :output)
      (setq 
       type (cadr args)
       args (cddr args)))
     ((vector list sequence)            ; others? eg string
      (setq 
       type (car args)
       args (cdr args)))
     ((arg :arg argument :argument parameter :parameter for :for)
      (setq    
       x (cadr args)
       args (cddr args)))
     (otherwise 
      (error 
       (tcl-or-elisp
        "Unknown ``keyword'' given to `image': ~s"
        "Unknown ``keyword'' given to `image': %s")     
       (car args)))))
  (unless
      (and D-supplied F-supplied)
    (error "must give both domain and formula to `image'"))
  (list 
   'map (list 'quote type) (list 'lambda (list x) F) D)))   
;
; Using mpacar when appropriate could be faster, since it doesn't have to check
; for type of argument, but just declaring type to compiler is usual way to get
; this added efficiency. 
;
; Possible extension: allow F to be a true function like `tanh' rather than an
; expression which gets construed as a function of its argument.
; To do this we might make it test for whether F is a 
; function or just an expression and act accordingly (but unfortuntely elisp
; lacks `functionp'.   But maybe don't need this extension very badly since
; mapcar itself is perspicuous in that instance.)
; An alternative is to use different keywords, with `function' @ a true fcn,
; and  `of' or `form' for usual case (for this reason, we do not have any such
; "function" keywords above).
;
; The idea of allowing use of just (image F D) with no keywords was dropped
; to avoid the danger of not catching errors of the sort,  (image of (* $ $)) ,
; i.e. an error of forgetting to put in the other clause.
;
; We can't use "?" to name the dummy variable since it has special meaning to
; elisp; we used `$' instead. 
  
;:----------------------------------------------------------
;: The macro kwa
;;--------------------------------
 
(localize-using-name "kwa" (ind b c to-key e H final)
(csynf   
(defmacro kwa (ind b c to-key e &rest H)
  "\
 A macro for simple looping, used in one of these three ways:
     (kwa j from 3 to   7  FORMS)    (j will be 3 4 5 6 7)
     (kwa j from 3 upto 7  FORMS)    (j will be 3 4 5 6)
     (kwa j from 3 while TEST FORMS)
 See the macro itself for some warning comments. \
  "
  (declare (ignore b))                 ; to stop compiler warning
  ;
  (case to-key
  ;-----------------------
  ;/case of inclusive `to'
  ;-----------------------
   ((to)			       
    `(csynf				; need csynf here to make it work!
       (varbind final ,e)
       (setq ,ind ,c)
       (while (<= ,ind final)
         ,@H
         (setq ,ind (1+ ,ind)))))
  ;-----------------------
  ;/case of exclusive `to'
  ;-----------------------
   ((upto)                             
    `(csynf
       (varbind final ,e)
       (setq ,ind ,c)
       (while (< ,ind final)
         ,@H
         (setq ,ind (1+ ,ind)))))
  ;----------------
  ;/case of `while'
  ;----------------
   ((while)
    `(csynf
       (setq ,ind ,c)
       (while ,e
         ,@H
         (setq ,ind (1+ ,ind))))) 
  ;---------------------------------------------
  ;/signal error if in none of the above 3 cases
  ;---------------------------------------------
   (t (error "Wrong syntax for `kwa': expected `to',`upto' or `while'"))))) )
 ;
 ; NOTES
 ;
 ; There is almost no need for this macro, since can use instead `loop',
 ; which exists in CL.  One possible reason to keep it, nevertheless, is that
 ; it is free of variables like `G6743'.
 ;
 ; The name is `kwa' rather than `for' because `edebug' didn't like latter.
 ;
 ; NB: crucial to use the INNER csynf's above for Clisp (not for elisp)
 ;
 ; This macro has exactly the problem described in the *info* documentation:
 ; the final value is evaluated each time thru, making it slower
 ;
 ; Improvement: extend to  (kwa j from a to b by c ....)
 ; where must handle the negative sign of c separately
 ; (can we make `by c' optional and not mess it all up? we could by changing to
 ; this format:  (kwa (j a b c) ....)
 ;
 ; Have localized `final' using `LUN', whence `(varbind final ,e)' is not
 ; really needed.  
 ;
 ; Have NOT localized the loop index since it is visible to the
 ; user.  However, one might want to do so using `(let ((j...' or by generating
 ; new symbols (on H as well))
 
;:: WARNING It takes a  L O N G  time to define `kwa' in clisp and elisp
 
;:-----------------------------------------------------
;: The arith-to-lisp family
;;----------------------------------------------------------------

;;; Next Still needs further checking and localization (see below)

(LUN "atl" (op L tail head)

(defmacro arith-to-lisp (&rest L) 
  " 
 Converts an expression in which the arithmetic operations `+' `-' `*' `/'
 have their usual syntax to one in which they have the lisp syntax.  It is not
 foolproof, perhaps; in particular you might get trouble if you use `+' etc. as
 variables, or as arguments to a macro situated inside this one.
  "
  (artlint L))

(defalias 'atl 'arith-to-lisp)

(defun artlint (L) 
  " 
 Merely the name for the internal function which does the real work of the
 macro arith-to-lisp.
  "
  (declare (special L head tail))
  ;
  (cond
   ;
   ((atom L) L)                         ; L is not a cons (so may be nil)
   ;
   ((null (cdr L))                      ; so L is a singleton
    (artlint (car L)))
   ;
   ((eq (car L) '+) (artlint (cdr L)))  ; begins with `+'
   ;
   ((eq (car L) '-)                     ; begins with `-'
    (cond
      ((= 2 (length L)) (list '- (artlint (cadr L))))              
      ((< 2 (length L)) (artlint (cons (list '-  (cadr L)) (cddr L))))
      (t (error "logically this message can never occur!"))))
   ;
   (t (let ((head nil) (tail nil))      ; list not beginning with + or -
        (declare (special head tail))     
        (cond 
         ((atl-locate '+) (list '+ (artlint head) (artlint tail)))
         ((atl-locate '-) (list '+ (artlint head) (artlint (cons '- tail))))
         ((atl-locate '*) (list '* (artlint head) (artlint tail)))
         ((atl-locate '/) (list '/ (artlint head) (artlint tail)))
          ;; (t (error " encountered an embedded non-arithmetic expression"))
           ; last removed to allow occurences inside of things like (length X)
         (t L))))))

(defun atl-locate (op) " An auxiliary function for arith-to-lisp." 
  (declare (special L head tail))  
  (setq tail (memq op L))
  (cond
   ((null tail) nil)
   (t
    (setq 
       head (subseq L 0 (- (length L) (length tail)))
       tail (cdr tail))
    (if 
        (or (null tail) (null head)) 
        (error
         "syntax error: invalid arithmetic operator at beginning or end"))
    op)))    
 )

 ; NOTES 
 ;
 ; 0. an atom is *any* non-cons-cell, including nil.
 ;
 ; 1. `atl-locate' returns nil if op not found, else op itself
 ;
 ; In elisp the `special' declarations only serve to stop
 ; compiler warnings.  In TCL perhaps they make the variable visible outside
 ; the function in which it is used ("indefinite extent")
 ;
 ; 2. TCL note.
 ;    Observe the use of (declare (special...)) both after the heading
 ;    of artlint AND after the (let ...).  It is needed in both places since
 ;    the 
 ;    `let' will otherwise create a local lexical binding for its arguments.
 ;    Is there some way to make the top declaration reach inside the `let'?

 ; HOW IT WORKS
 ;
 ;   At all stages the working expression is a list of elements to be thought
 ;   of as "operators" and "operands".  (But the argument to the macro
 ;   itself is just just an arithmetic expression not requiring overall
 ;   enclosure parentheses). 

 ; MOTIVATION
 ;
 ; We want to be able at least to evaluate things like
 ;
 ;    - A + 1/2 (B - C) / (d % x)
 ;
 ; But maybe not wise to do this, since something like
 ; 
 ;            (atl 4 + (length L))
 ; 
 ; might confuse it to think `length'is a vble, not a function.  Here there is
 ; actually no ambiguity, can there ever be?  Probably not, unless you decide 
 ; to use `+' etc. as a variable, or as an argument to a macro.  So don't be
 ; deterred by this worry.
 ;

;;; NEEDED:
 ;
 ;; Should do more checking for invalid syntax. (in particular could convert
 ;; the commented (;;) out error message into a warning.
 ;
 ;; Should sort out exactly what expressions are handled by atl and which are
 ;; not (see tests for some examples), specifically when functions are
 ;; evaluated within a call to atl. 
 ;
 ;; Tests are not extensive yet.

;; POSSIBLE IMPROVEMENTS
;;
;; 1. make (a - b) become (- a b) rather than (+ a (- b)).  This can be done by
;;    handling the minus signs from right to left rather than left to right.
;;    Then it might also be unnecessary to treat the case of length=2
;;    separately near the beginning of artlint.
;;
;; 2. make (a + b + c) become (+ a b c) rather than (+ a (+ b c))
;;
;; 3. could even combine these to make (a - b - c + d - e + f) become
;;    (- (+ a b f) (+ b c e))
;;
;; 4. Handle more infix operators:  ^ == **  and % and maybe ! for "choose"
;;    also comparison operators like > < = >= etc.
;;
;; 5. The auxiliary functions could be written inside the main one using
;;    `flet'or friends. On the other hand, we
;;    don't want them to be defined each time the macro is called!!
 
;:------------------------------------------------------ 
;: Other macros 

(when *elisp*
(defmacro retrieve-error-messages-from (&rest body)
  " 
 Will return the error information (if an error occurs), allowing 
 retrieval of error messages which overrun the echo area.  Actually this is no
 longer needed because in emacs 19.30 and later, the messages are saved in the
 buffer *Messages*
  " 
  `(condition-case msg
       (progn ,@ body)
        (error msg))) )
  ;  
  ; Once upon a time, used this macro to see error messages which get lost
  ; because the echo area is too small.  
  ; It doesn't retrieve the error message as such, but the symbol naming the
  ; error, together with the other data which would appear in the message.
  ;
  ; equivalent code not using backquote:
  ;
  ;   (list 'condition-case 'msg
  ;      (cons 'progn body)
  ;      (list 'error 'msg)))



;; COMMENTS ON TIMER FAMILY
;; The time macros are in a mess as of now.
;; putting more or less (progn...) in things changes their behavior!!!

(when *elisp*
  (localize-using-name "timer-family" (body form t_0 t_1 time-triple z)
  (defvar t_0)
  (defvar t_1)
 ;
(defun get-time ()  
  "\
 Returns current time in seconds from some meaningless reference time.
 (Must rewrite this if elapsed times will be longer than 99999 seconds,
 which is around 12 days.)"  
  (let
      ((time-triple (current-time)))
       (+
        (cadr time-triple)  
        (/ (caddr time-triple) 1000000.0))))
        ;
(defun pp-elapsed-time (z) 
  "\
 Argument should be in seconds.  constructs a sentence saying that this is the
 amount of time that elapsed and calls princ on it (more or less)
  "
  (cond
   ((< z 0) (message "Negative elapsed time of %f seconds!" z))
   (t
    (princ
     (format "%0.3g sec elapsed (NOT cpu time)" z)) )))  
     ; 
(defun time-evaluation (form) 
  "\
 A very crude timer that only measures elapsed real time, not CPU time. 
 Evals its argument and returns the elapsed time in seconds. 
 BEWARE: on bananoid this behaves differently when enclosed in `progn'!
  "
  (setq t_0 (get-time))
  (eval form)
  (setq t_1 (get-time))
  (- t_1 t_0))  
   ;
   ;  The use of `progn' somehow stabilizes the timing (on bananoid)
   ;
(defmacro time (body) 
  "\
 A very crude timer, which only measures elapsed real time.
 Prints a message with elapsed time.   (Actual return value is `t')"
  (pp-elapsed-time (time-evaluation body))
  t)  
       ;
(defmacro Time (&rest body)
  " First collect garbage then time (still only elapsed real time)"
  `(progn
     (garbage-collect)
     (time
      ,@body)))          ))
 
;:----------------------------------------------------------
;: Macros specifically for TCL (for compatibility with elisp) 
;------------------
(when *elisp* (defalias 'elisp-if 'if))

(when *gcl*

  (defmacro elisp-if (A B &rest C)
  " This is the elisp form of `if', for use with old elisp code.  It allows
 multiple forms in the ``else'' clause."
  (list 'if A B (cons 'progn C)))

  (defmacro while (test &rest body) " Reproduces this macro for sun lisp"
    (list 'loop 'while test 'do (cons 'progn body)))

  ; (if *clisp*
  ; (defmacro while (test &rest forms) " Reproduces this macro for clisp"
  ;   `(do () ((not ,test)) ,@forms)))

  )
;:----------------------------------------------------------
