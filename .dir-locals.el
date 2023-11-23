((nil .
      ((eval .
             (progn
               (defun fql/tmux-pane-cmd (pane cmd)
                 "Runs the given shell command in a subshell inside a tmux pane."
                 (interactive)
                 (let* ((resolved-pane (concat "fastqtl-lmm:" pane))
                        (resolved-cmd (format "'%s'" cmd))
                        (cmd-parts (list "tmux"
                                         "send-keys"
                                         "-t"
                                         resolved-pane
                                         resolved-cmd
                                         "C-m")))
                   (shell-command (format "tmux clear-history -t %s" resolved-pane))
                   (shell-command (mapconcat 'identity cmd-parts " "))))

               (defun get-file-in-project (filename)
                 "Gets absolute path to file inside the project"
                 (interactive "P")
                 (expand-file-name filename (projectile-project-root)))

               (defun fql/run-sh-scratch ()
                 "Runs the shell scratch file"
                 (interactive)
                 (fql/tmux-pane-cmd "0.0" "(clear; ./tmp/scratch.sh)"))

               (defun fql/kill-dev-env ()
                 "Kills all processes and tmux"
                 (interactive)
                 (fql/tmux-pane-cmd "0.0" "C-c")
                 (fql/tmux-pane-cmd "0.0" "tmux kill-session"))

               (global-set-key (kbd "<f1>") 'fql/run-sh-scratch)
               (global-set-key (kbd "<f12>") 'fql/kill-dev-env)

               )))))
