# Release process

MetaQuest publishes Linux-supported Python distributions through GitHub Actions
and PyPI Trusted Publishing. No long-lived PyPI API token is stored in GitHub.

## One-time repository configuration

1. In GitHub, create environments named `testpypi` and `pypi`.
2. Require a maintainer approval for both environments. Restrict `pypi` to
   protected tags matching `v*` if the repository plan exposes that control.
3. Protect `main` and require these CI checks before merging:
   `Python 3.10`, `Python 3.11`, `Python 3.12`,
   `Direct dependency licenses`, `Build distributions`, and
   `Conda integration (Python 3.11)`.
4. On TestPyPI, configure a trusted publisher for:
   - owner: `dpatel511`
   - repository: `metaquest`
   - workflow: `testpypi.yml`
   - environment: `testpypi`
5. On PyPI, configure a trusted publisher for the same owner and repository,
   using workflow `release.yml` and environment `pypi`. A pending publisher can
   reserve the project name for its first upload.

## TestPyPI rehearsal

Run the **Publish to TestPyPI** workflow manually. Approve the `testpypi`
environment deployment, then inspect both the project page and uploaded wheel
and source distribution. Versions on package indexes are immutable, so bump the
pre-release serial before repeating an upload.

## Production release

1. Update `src/metaquest/config.py` to the canonical PEP 440 version and add the
   corresponding changelog entry.
2. Merge the release preparation PR after all required checks pass.
3. Create a signed or annotated tag with the same canonical version prefixed by
   `v`, for example `v2.0.0a1`.
4. Create and publish a GitHub Release from that tag.
5. Approve the `pypi` environment deployment.

Publishing the GitHub Release triggers `release.yml`. The workflow verifies
that the tag and package versions match, builds and checks the wheel and source
distribution in one job, and gives only the final upload job an OIDC token.
Draft releases do not publish to PyPI.
