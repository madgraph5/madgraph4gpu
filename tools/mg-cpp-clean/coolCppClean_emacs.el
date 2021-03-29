;;============================================================================
;;  LHCb customisations
;;============================================================================

;; Avoid double indentation at statement opening
;; NB This has no effect here, so it is added to the emacs-cppclean function below
;; NB Alternatively put it in ~/.emacs and load it via 'emacs --batch -u $USER...'
;;(c-set-offset 'substatement-open 0)  

;;(setq c++-fill-column 130 )        ;; not needed?
;;(setq comment-line-start "//")     ;; not needed?
(setq-default indent-tabs-mode nil)  ;; Never insert an ASCII TAB
(setq tab-width 2)                   ;; TAB is two spaces

;; Deletes the whitespace at the end of any line
(defun delete-eol-whitespace ()
  "Deletes the whitespace at the end of any line"
  (interactive)
  (save-excursion
    (goto-char (point-min))
    (while (< (point) (point-max))
      (forward-line 1)
      (end-of-line)
      (delete-horizontal-space)
      )
    )
  )

;;============================================================================
;;  COOL customisations
;;============================================================================

(defun myindent-buffer ()
  "Indent buffer excluding comments."
  (interactive)
  (save-excursion
    (let (pmin pmax commentstart commentend skipstart skipend pskip indcstart pindc))
    (goto-char (point-min))
    (setq pmin (point-min))
    (setq pmax (point-min)) 
    (setq commentstart "/*")
    (setq commentend "*/")
    (setq skipstart (concat commentstart " COOLCPPCLEAN-NOINDENT-START " commentend))
    (setq skipend (concat commentstart " COOLCPPCLEAN-NOINDENT-END " commentend))
    (setq pskip nil)
    (setq indcstart (concat commentstart "*")) ;; INDENT COMMENTS STARTING WITH "/**"!
    (message "Indent buffer (excluding comments) from %d to %d..." (point-min) (point-max))
    (while (< pmax (point-max))
      ;; 1. LOOK FOR START OF COMMENT (OR START OF NOINDENT SECTION)
      (setq pindc nil)
      (setq pmax (search-forward commentstart (point-max) t))
      (if pmax 
          ;; Start of comment found: look for indented comments and noindent sections
          (progn
            ;;(message "DEBUG: '%s' found at %d" commentstart pmax)
            ;; Start of comment found: is this the start of an indented comment?
            (backward-char 2)
            (setq pindc (search-forward indcstart (point-max) t) )
            ;;(if pindc (message "DEBUG: '%s' found at %d" indcstart pindc))
            (if pindc
                ;; The start of an indented comment was found: is it here?
                (if (not (= (- pindc pmax) (- (length indcstart) 2))) 
                    (progn 
                      (setq pindc nil)
                      ;;(message "DEBUG: '%s' ignored" indcstart)
                      (goto-char pmax) )
                  ) 
              )
            ;;(if pindc (message "DEBUG: '%s' confirmed at %d" indcstart pindc))
            ;; Start of comment found: is this the start of a noindent section?
            (if pindc ()
              (progn
                (backward-char 2)
                (setq pskip (search-forward skipstart (point-max) t) )
                ;;(if pskip (message "DEBUG: '%s' found at %d" skipstart pskip)) 
                )
              (if pskip 
                  ;; The start of a noindent section was found: is it here?
                  (if (not (= (- pskip pmax) (- (length skipstart) 2))) 
                      (progn 
                        (setq pskip nil)
                        ;;(message "DEBUG: '%s' ignored" skipstart)
                        (goto-char pmax) )
                    ) 
                )
              ;;(if pskip (message "DEBUG: '%s' confirmed at %d" skipstart pskip)) 
              )
            ) 
        ;; No comment found: break the loop
        (progn 
          ;;(message "DEBUG: '%s' not found: break the loop" commentstart)
          (setq pmax (point-max)) 
          (setq pskip nil) ) )
      ;; Indent until the start of the comment (or the end of the buffer)
      ;;(if (not pindc) (message "Indent from %d to %d" pmin pmax))
      (if (not pindc) (indent-region pmin pmax nil) )
      ;; 2. LOOK FOR END OF COMMENT (OR END OF NOINDENT SECTION)
      (if (or pindc (>= pmax (point-max)))
          ;; No start of comment/noindent, or indented comment: no need to match end
          (progn 
            ;;(message "DEBUG: no need to find matching '%s'" commentend) 
            )
        ;; Start of comment or noindent found: look for the matching end
        (progn
          (let (tagend))
          ;; What is the matching end: are we inside a comment or a noindent section?
          (setq tagend (if pskip skipend commentend))
          (setq pmin (search-forward tagend (point-max) t))
          (if (not pmin) 
              ;; The end was not found (BAD C++ CODE!): break the loop
              (progn
                (setq pmax (point-max))
                (message "WARNING: matching '%s' not found" tagend) )
            ;; The end was found
            ;;(message "DEBUG: matching '%s' found at %d" tagend pmin) 
            )
          )
        )
      )
    (message "Indent buffer (excluding comments) from %d to %d... done" (point-min) (point-max))
    )
  )

;;============================================================================

(defun clean ()
  "Modified LHCb clean() function, avoiding the indentation of comments."
  (interactive)
  (save-excursion
    (delete-eol-whitespace)
    ;; TAB replaced by space
    (goto-char (point-min))
    (while (and (re-search-forward "\t" nil t ) 
		(< (point) (point-max)))
      (replace-match " " nil nil) )
    ;; CR replaced by nothing (but "\r" should be C-m making this a binary file)
    ;;(goto-char (point-min))
    ;;(while (search-forward (string "\r") nil t)
    ;;		(< (point) (point-max)))
    ;;(replace-match "" nil nil) )
    ;; Alternative dos2unix
    (set-buffer-file-coding-system 'utf-8-unix)
    ;; Re-indent the buffer
    (myindent-buffer)
    )
  )

;;============================================================================

(defun emacs-cppclean ()
  "Clean C++ code."
  (widen)
  (message "Clean C++ code from %d to %d..." (point-min) (point-max))
  (c++-mode)                          ;; NB (c++-mode) must be BEFORE (c-set-offset..)
  (c-set-offset 'substatement-open 0) ;; Avoid double indentation at statement opening
  (clean)
  (save-buffer)
  (message "Clean C++ code from %d to %d... done" (point-min) (point-max))
  )

;;============================================================================

(defun emacs-pyclean ()
  "Clean Python code."
  (widen)
  (message "Clean Python code from %d to %d..." (point-min) (point-max))
  (python-mode) ;; NB better call (python-mode) BEFORE (untabify...)
  (untabify (point-min) (point-max))
  (clean)
  (save-buffer)
  (message "Clean Python code from %d to %d... done" (point-min) (point-max))
  )

;;============================================================================

;; Disable Python smart indentation and hardcode offset to 4:
;; by default smart indentation was on and guessed offsets between 2 and 8,
;; so that if an existing file had offset 8, M-x clean did not have any effect.
;; Use M-x describe-variable py-indent-offset to see the buffer local value...
;; See http://svn.python.org/projects/python/trunk/Misc/python-mode.el
;; See http://jesselegg.com/archives/2010/02/25/emacs-python-programmers-part-1
(setq-default py-indent-offset 4)
(setq-default py-smart-indentation nil)

;; Just to be sure, disable indent-tabs-mode whenever in python mode,
;; otherwise some tabs may be left in python code instead of spaces.
;; See http://stackoverflow.com/questions/4251159/set-python-indent-to-2-spaces-in-emacs-23
;; PS: This was DISABLED - instead, call (untabify) before (clean)...
;; Otherwise you may get 'Bad indentation errors' from emacs -batch!
;;(add-hook 'python-mode-hook
;;              (lambda ()
;;                (setq indent-tabs-mode nil)
;;                (setq tab-width 4)))

;;============================================================================

