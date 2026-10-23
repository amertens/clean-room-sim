# cleanTMLE 0.5.0: release notes and steps (not yet published)

Prepared 2026-10-20. Nothing here has been pushed, tagged or released. Publishing needs Andrew's word.

## To publish

From `cleanTMLE/`, in order:

```bash
git push origin main
```

```bash
git tag -a v0.5.0 -m "cleanTMLE 0.5.0"
```

```bash
git push origin v0.5.0
```

```bash
gh release create v0.5.0 --title "cleanTMLE 0.5.0" --notes-file ../docs/cleantmle-0.5.0-release-notes-body.md
```

## Release notes (the body of the GitHub release)

The text is in `cleantmle-0.5.0-release-notes-body.md`, next to this file.
