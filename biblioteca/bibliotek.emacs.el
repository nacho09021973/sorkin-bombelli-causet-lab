;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;                                                          ;;;;;;;;
;;;;;;  All files in this directory are copyright 1997, 1998,   ;;;;;;;;
;;;;;;  1999 by Rafael D. Sorkin.  All rights reserved.         ;;;;;;;;
;;;;;;                                                          ;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                            ;Time-stamp:<Nov 05 1998 03:38:34 (13889 25610)>

;=========================================================================
; Functions to enhance the operation of emacs (as an editor)
; There are also a few strays here like `macroexpand-region'
;
;   They could go into { .emacs } but are put here FOR AUTOLOADING
;   instead (because they use csynf, for example)
;
;   WE DON'T COMPILE THEM (neither are the functions in .emacs compiled).
;=========================================================================
; Index of Functions and vbles 
; 
;   utupe
;   tunza (MACHINE DEPENDENT)
;   lisp-insert-comment
;   macroexpand-region
;
;   andika-file-kwa-rcp
;   nakili-file-umojani (andika-file-umojani is commented out)
;   andika-file-sunix
;   
;-----------------------------------------------------------------------

(require 'preparations "~/lisp/preparations.el")

(deff utupe (&optional overwrite)
  "\
 This is intended to provide `utupe' from within dired.  It will move the
 file on the current line to the directory ~/.Taka, but will not overwrite 
 a file already there, unless the optional arg `overwrite' is nonnil.  \
  "
  (interactive "P")
  (require 'preparations "~/lisp/preparations.el")
  (varbind
    x (dired-get-filename)
    y (concat "~/.Taka/" (dired-get-filename 'no-dir))
    ow (o not not overwrite))
  (when
    (y-or-n-p 
     (format " ? %s ?   " (dired-get-filename 'no-dir)))
    (dired-rename-file x y ow)))
 ;
 ; Notes
 ;
 ;  The third argument to `dired-rename-file' is a flag, whose meaning is not
 ;  documented.  Empirically, if it is nil you won't overwrite another file of
 ;  the same name, if it is t then you will.  Are there other values?  
 ;
 ;  Instead of providing an optional argument, as we have done, you
 ;  could get it to ASK whether to overwrite, by use of error handling.
 ;  Or (more simply?) could just check for the file there, and ask first
 ;  whether to overwrite if it's present.
 ;
 ;  It is important in `interactive' that you use big "P" because using little
 ;  "p" makes a void argument turn into the number 1 (!)
 ;
 ;  The alternative approach to use `dired-do-shell-command' is bad, because it
 ;  won't update the file list.

;---------------------------------------------------------------

;;; should we make this uuencode the file first?

(deff tunza ()
  " Yatunza the current buffer - SUHEP version."
  (interactive)
 ;--------------------------------------------------------
 ;/specify where to send the backup copies (SUHEP VERSION)
 ;--------------------------------------------------------
  (varbind 
   destinations
   (list
    "sorkin@vega.npac.syr.edu"
    "sorkin@nuclecu.unam.mx"
    ; "fhusein@sirius.fis.cinvestav.mx"
    ; "sorkin@suhep.phy.syr.edu"
    ; "rdsorkin@mailbox.syr.edu"
    ))
 ;-----------------------------------------------------------------------
 ;/specify buffer to take output from shell command and name of temp file
 ;-----------------------------------------------------------------------
  (varbind 
   command-out "*kutunza*"
   tempfile "*kutunza*.tmp")
 ;-----------------------------------------------------
 ;/get name of file and extract its last two components
 ;-----------------------------------------------------
  (varbind 
   full-filename buffer-file-name
   short-name
   (subseq 
    full-filename
    (1+ (string-match "/[^/]*/[^/]*$" full-filename))))
 ;-------------------
 ;/update time-stamp 
 ;-------------------
  (time-stamp)
 ;----------------------------------------
 ;/write buffer contents to temporary file
 ;----------------------------------------
  (write-region (point-min) (point-max) tempfile)
 ;---------------------------
 ;/ask whether to save buffer
 ;---------------------------
  (if (y-or-n-p "Save file locally? ") (save-buffer))
 ;----------------------------------
 ;/Mail copies out to `destinations'
 ;----------------------------------
  (loop 
   for kuda in destinations do
   (shell-command
    (concat
     "< " tempfile
     " Mail "
     " -s " (concat "'" "[backup]" short-name "'")
					; " -s " (concat "[backup]" short-name)
     " " kuda) 
    command-out))
 ;----------------------
 ;/delete temporary file
 ;----------------------
  (delete-file tempfile)
 ;-------------------------------------------------------------------
 ;/finally issue completion message, and show output if there was any
 ;-------------------------------------------------------------------
  (message "imetunzwa: %s " short-name)
  (if (get-buffer command-out) 
      (switch-to-buffer (get-buffer command-out))))
  ;
  ; NOTE
  ; 
  ; It's not clear how to get the return status of the command (it seems NOT to
  ; be the value returned by `shell-command', which rather seems to be `t' no
  ; matter what happens).  So we can't check the return status directly to
  ; detect problems.  
  ; Output from the command should appear in buffer "*Shell Command Output*".
  ; If an error message comes from the shell it should appear in 
  ; buffer "*kutunza*" and be displayed.


(deff lisp-insert-comment () 
  " insert comment line into lisp code before current line"
  (interactive)
 ;----------------------------------------------------------
 ;\get the text of the comment and make a border line for it
 ;----------------------------------------------------------
  (varbind 
   comment-line (read-from-minibuffer "comment: " ";/")
   border-line (make-string (1- (length comment-line)) ?-))   
 ;------------------
 ;/insert top border
 ;------------------
  (beginning-of-line)
  (open-line 1)
  (indent-for-tab-command)
  (backward-delete-char-untabify 1)
  (insert (concat ";" border-line))
  (newline)
 ;----------------------
 ;/insert comment itself
 ;----------------------
  (indent-for-tab-command)
  (backward-delete-char-untabify 1)
  (insert comment-line)
  (newline)
 ;---------------------
 ;/insert bottom border
 ;---------------------
  (indent-for-tab-command)
  (backward-delete-char-untabify 1)
  (insert (concat ";" border-line)))

(deff macroexpand-region ()
  " REPLACE region by its pretty-printed macro expansion "
  (interactive)
  (varbind H (read (buffer-substring (region-beginning) (region-end))))
  (setf 
   (buffer-substring (region-beginning) (region-end))
   (pp (macroexpand H))))

(deff andika-file-kwa-rcp (destination)
  "\
 Use `rcp' to copy the file associated with the current buffer to 
 the machine specified.
 The argument should look like \"rdsorkin@umoja.syr.edu\"
 BEWARE: will overwrite the remote file, if there is one.
  "
  (varbind filename (o abbreviate-file-name buffer-file-name))
  (shell-command
   (format
    "rcp %s %s:%s &" filename destination filename)
   "*Remote Copy Output*")
  (display-buffer "*Messages*")
  (switch-to-buffer "*Remote Copy Output*"))
  ;; seems impossible to get top line of buffer to show on display!


(deff nakili-file-umojani ()
  "\
 Copy buffer (destructively) to file of same name on umoja disks.
 Use this only when you are logged into umoja ``as if on suhep''.  
 If you are really on suhep this won't work!
  "
  (interactive)
  ;
  (varbind
   original-home "/home1/sorkin/"
   copy-home "/users/rdsorkin/")
  ;
  (unless (equal system-name "umoja.syr.edu")
    (error "Can only use this from umoja"))
  (unless
      (equal (subseq buffer-file-name 0 14) original-home)
    (error "Only use this if you are ``pretending to be at suhep''"))
 ;----------------------------------
 ;/get name of file @ current buffer and form further names from it
 ;----------------------------------
  (varbind 
    filename (subseq (abbreviate-file-name buffer-file-name) 2)
    tempname (format "%s%s.tmp" copy-home filename)
    copy-name (format "%s%s"    copy-home filename))
 ;----------------------------------------------
 ;/protect from overwriting eigen-files of umoja
 ;----------------------------------------------
  (if (equal filename ".bashrc.local")
      (error "Don't do this with { %s } !" filename))
 ;------------------------------------------------------------
 ;/write to temp file in counterpart directory on umoja disks
 ;------------------------------------------------------------
  (write-region (point-min) (point-max) tempname)
 ;---------------------------------
 ;/rename temp file to correct name
 ;---------------------------------
  (rename-file tempname copy-name 'overwrite)
  (message "written: %s" copy-name))
 ;
 ; NOTE
 ; Reason for copying first and then moving copy is to protect previous
 ; version from being erased.  If a writing error occurs due to full disks,
 ; then this should only affect the temp file: execution should stop at that
 ; point and old file will not be touched.

  ;(defun andika-file-umojani () 
  ;  " Copy this buffer's file to umoja (destructively)"
  ;  (interactive)
  ;  (andika-file-kwa-rcp "rdsorkin@umoja.syr.edu"))
  ;
  ; doesn't work because `rsh' is confused by output from .bashrc et al (this
  ; is a known problem mentioned in the man page)

(defun andika-file-sunix () 
  " Copy this buffer's file to gamera (destructively)"
  (interactive)
  (andika-file-kwa-rcp "rdsorkin@gamera.syr.edu"))

