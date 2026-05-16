# Project Guidance

This repository is an archive of Rafael D. Sorkin materials, including Emacs Lisp source under `biblioteca/`, archived papers in `papers/`, and PostScript figures in `PostScript/`.

## Working Rules

- Prefer reading the Lisp sources and inline documentation before making changes.
- Treat PDF, PS, and compiled artifacts as reference material unless the task explicitly requires editing them.
- Keep edits minimal and localized to the relevant source file.
- Preserve the existing file naming scheme for historical materials.
- Do not rewrite or normalize the archived documentation style unless asked.

## Useful Entry Points

- `biblioteca/README` for the package overview and installation notes.
- `biblioteca/Some.more.documentation` for additional background.
- `biblioteca/preparations.el` and `biblioteca/preparations.gcl` for load/setup flow.

## Editing Workflow

- Start with the source file that owns the behavior:
  - `biblioteca/bibliotek.poset.el` for poset operations.
  - `biblioteca/bibliotek.general.el` for general list/set helpers.
  - `biblioteca/bibliotek.float.el` for numeric helpers.
  - `biblioteca/bibliotek.macros.el` for macro and syntax support.
- Search for the function name before editing so related helpers and docstrings stay consistent.
- If a behavior change affects multiple files, update the lowest-level helper first, then the callers.
- Keep public function names and argument order stable unless the task explicitly asks for an API change.
- Preserve the elisp/Common Lisp dual-use intent when touching shared code.

## Validation

- Prefer small, local sanity checks over broad rewrites.
- When possible, inspect the surrounding docstrings and examples to confirm the intended semantics.
- If the change is documentation-only, avoid code edits and keep the wording close to the existing style.

## Useful Commands

- `rg -n "pattern" biblioteca/` to find definitions, docstrings, and cross-references.
- `rg --files biblioteca/` to list the source and documentation files in the package.
- `sed -n '1,220p' biblioteca/<file>` to inspect the start of a Lisp source file or document.
- `find biblioteca -maxdepth 1 -type f | sort` to review the package layout quickly.
- `file biblioteca/<file>` when you need to confirm whether something is source, text, PDF, or compiled output.

## File Map

- `biblioteca/bibliotek.poset.el`: core poset, order, and relation operations.
- `biblioteca/bibliotek.extras.el`: additional poset helpers and higher-level utilities.
- `biblioteca/bibliotek.general.el`: list, set, and general-purpose helpers.
- `biblioteca/bibliotek.float.el`: numeric and floating-point helpers.
- `biblioteca/bibliotek.macros.el`: macros and language-level conveniences.
- `biblioteca/bibliotek.constants.el`: shared constants and configuration values.
- `biblioteca/bibliotek.emacs.el`: Emacs-specific integration.
- `biblioteca/bibliotek.elisp.el` and `biblioteca/bibliotek.elisp.patch.el`: Elisp-oriented compatibility material.
- `biblioteca/bibliotek.TCL.gcl`: Common Lisp / GCL-oriented support.

## Validation Notes

- If Emacs is available, a narrow `--batch` load of `biblioteca/preparations.el` is a reasonable smoke test for Lisp edits.
- Prefer validating the smallest affected path first, then expand only if the change crosses module boundaries.
- When no executable validation is practical, at least confirm the affected docstrings, examples, and call sites stay aligned.

## Suggested Smoke Tests

- `emacs --batch -l biblioteca/preparations.el` to load the package setup non-interactively.
- `emacs --batch --eval '(load-file "biblioteca/preparations.el")'` if you need an explicit load path.
- `emacs --batch --eval '(progn (load-file "biblioteca/preparations.el") (message "loaded"))'` for a minimal success check.

## Search Patterns

- `rg -n "^\\s*\\(defun\\|defmacro\\|deff\\)" biblioteca/bibliotek.poset.el` to list function-like definitions.
- `rg -n "^[[:space:]]*past\\b|^[[:space:]]*future\\b|^[[:space:]]*interval\\b" biblioteca/` to locate core poset APIs and references.
- `rg -n "count-.*for|t-close|sortedp|level\\b" biblioteca/` to find related helpers by naming convention.
- `rg -n "describe-function|docstring|TODO|FIXME" biblioteca/` to review embedded documentation and maintenance notes.

## Typical Changes

- For a new or revised poset operation, start in `biblioteca/bibliotek.poset.el`, then update any related notes in `biblioteca/Some.more.documentation` if needed.
- For helper behavior shared across modules, edit `biblioteca/bibliotek.general.el` or `biblioteca/bibliotek.float.el` before changing callers.
- For macro or syntax changes, edit `biblioteca/bibliotek.macros.el` and verify the resulting forms still work in both elisp and Common Lisp contexts.
- For Emacs-specific integration, keep the change in `biblioteca/bibliotek.emacs.el` so the core libraries remain portable.
- For documentation-only fixes, prefer `biblioteca/README`, `biblioteca/Some.more.documentation`, or the relevant file header instead of code files.

## Conventions

- Use ASCII when adding new text unless an existing file already uses non-ASCII characters.
- If you modify Lisp code, check nearby docstrings and comments for the intended behavior.
- Avoid broad refactors; this codebase appears to be historical and documentation-heavy.
