
# Contributing to V-pipe

A big welcome and thank you for considering contributing to V-pipe! It’s people like you that make it a reality for users in our community.

Reading and following these guidelines will help us make the contribution process easy and effective for everyone involved. It also communicates that you agree to respect the time of the developers managing and developing these open source projects. In return, we will reciprocate that respect by addressing your issue, assessing changes, and helping you finalize your pull requests.

## Quicklinks

* [Getting Started](#getting-started)
    * [Issues](#issues)
    * [Pull Requests](#pull-requests)
* [Getting Help](#getting-help)

## Getting Started

Contributions are made to this repo via Issues and Pull Requests (PRs). A few general guidelines that cover both:

- Search for existing Issues and PRs before creating your own.
- We work hard to makes sure issues are handled in a timely manner but, depending on the impact, it could take a while to investigate the root cause. A friendly ping in the comment thread to the submitter or a contributor can help draw attention if your issue is blocking.

### Issues

Issues should be used to report problems with the V-pipe workflow, request a new feature, or to discuss potential changes before a PR is created. When you create a new Issue, a template will be loaded that will guide you through collecting and providing the information we need to investigate.

If you find an Issue that addresses the problem you're having, please add your own reproduction information to the existing issue rather than creating a new one. Adding a [reaction](https://github.blog/2016-03-10-add-reactions-to-pull-requests-issues-and-comments/) can also help be indicating to our maintainers that a particular problem is affecting more than just the reporter.

### Pull Requests

PRs to our workflow are always welcome and can be a quick way to get your fix or improvement slated for the next release. In general, PRs should:

- Target our staging branch: [rubicon](https://github.com/cbg-ethz/V-pipe/tree/rubicon)
- Only fix/add the functionality in question **OR** address wide-spread whitespace/style issues, not both.
- Add unit or integration tests for fixed or changed functionality (if a test suite already exists).
  - Or at least provide a minimalist example dataset
- Address a single concern in the least number of changed lines as possible.
- Include documentation in the repo or on our `docs/` directory.

For changes that address core functionality or would require breaking changes (e.g. a major release), it's best to open an Issue to discuss your proposal first. This is not required but can save time creating and reviewing changes.

In general, we follow the ["fork-and-pull" Git workflow](https://github.com/susam/gitpr)

1. Fork the repository to your own Github account
2. Clone the project to your machine
3. Create a branch locally with a succinct but descriptive name
4. Commit changes to the branch
5. Following any formatting and testing guidelines specific to this repo
   - We rely on [snakefmt](https://github.com/snakemake/snakefmt) for Snakemake files
   - We use [Mega-Linter](https://megalinter.io) for the remaining files (Python (Black), Jupyter (Jupyfmt), Markdown (Markdownlint), Bash (Shellcheck), Perl (Perlcritic), Docker (Hadolint))
   - Ask us for help if you have trouble linting your code
6. Push changes to your fork
7. Open a PR in our repository and follow the PR template so that we can efficiently review the changes.

## Getting Help

Join us in the [V-pipe Gitter channel](https://gitter.im/V-pipe/community) (also [accessible over matrix](https://matrix.to/#/#V-pipe_community:gitter.im?utm_source=gitter) from your favorite client) and post your question there to reach out the devs.
For further inquiries, you can also contact the V-pipe Dev Team by opening a ticket at [v-pipe@bsse.ethz.ch](mailto:v-pipe@bsse.ethz.ch).
