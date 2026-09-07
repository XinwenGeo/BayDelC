# Release checklist

The package is prepared through the final local stage before GitHub upload.

- [x] Public API, model registry, S3 methods, and manual pages
- [x] Compact model assets generated without refitting
- [x] Model manifest and extraction/validation scripts
- [x] README, examples, tutorials, citation, license, and contribution policy
- [x] Cross-platform GitHub Actions R CMD check workflow
- [x] Local source build, clean installation, examples, and smoke tests
- [x] Numerical comparison against the archived full model bundles
- [x] Main/sensitivity selection and validation of every distributed asset
- [x] Package-native article Figure 2 and Figure 3 examples
- [x] Inert calibration reference scripts, Stan files, and source tables
- [ ] Replace the in-preparation article citation with final DOI/journal metadata
- [ ] Create the GitHub repository and push the reviewed package directory
- [ ] Confirm GitHub Actions passes on macOS, Windows, R release, and R devel
- [ ] Create a versioned GitHub release after publication

Do not add raw MCMC chain files, compiled Stan artifacts, submitted figure
outputs, or the full paper-specific plotting tree to this repository.
