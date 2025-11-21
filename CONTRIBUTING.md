# Contributing to Rare Variant GWAS Lollipop Plot Generator

Thank you for your interest in contributing! This document provides guidelines for contributing to this project.

## Ways to Contribute

- **Report bugs:** Open an issue describing the bug and how to reproduce it
- **Suggest enhancements:** Open an issue describing your enhancement idea
- **Submit code:** Fork the repository, make changes, and submit a pull request
- **Improve documentation:** Help us make the README and documentation clearer
- **Share examples:** Contribute example datasets or use cases

## Getting Started

1. **Fork the repository** on GitHub
2. **Clone your fork:**
   ```bash
   git clone https://github.com/YOUR_USERNAME/Rare-variant-GWAS-lollipop-plot.git
   cd Rare-variant-GWAS-lollipop-plot
   ```

3. **Create the development environment:**
   ```bash
   conda env create -f environment.yml
   conda activate lollipop
   ```

4. **Create a new branch:**
   ```bash
   git checkout -b feature/your-feature-name
   ```

## Code Guidelines

### R Code Style

- Follow standard R coding conventions
- Use meaningful variable names
- Comment complex logic
- Keep functions focused and single-purpose
- Use data.table or tidyverse consistently

### Testing Your Changes

Before submitting a pull request:

1. **Test basic functionality:**
   ```bash
   Rscript lollipop_maker.R --symbol GCDH
   ```

2. **Test with different parameters:**
   ```bash
   Rscript lollipop_maker.R --symbol GCDH --qvalue 1e-5 --output test_output/
   ```

3. **Verify outputs:**
   - Check that plots are generated correctly
   - Verify the TSV output has the expected format
   - Ensure no errors or warnings appear (unless expected)

### Documentation

- Update README.md if you change functionality
- Add comments to complex code sections
- Update QUICK_REFERENCE.md if you add new features
- Include usage examples for new features

## Submitting Changes

1. **Commit your changes:**
   ```bash
   git add .
   git commit -m "Brief description of your changes"
   ```

2. **Push to your fork:**
   ```bash
   git push origin feature/your-feature-name
   ```

3. **Open a Pull Request:**
   - Go to the original repository on GitHub
   - Click "New Pull Request"
   - Select your fork and branch
   - Describe your changes clearly
   - Reference any related issues

## Pull Request Guidelines

- **Clear description:** Explain what your PR does and why
- **One feature per PR:** Keep PRs focused on a single feature or fix
- **Test your code:** Ensure it works with the example data
- **Follow the code style:** Match the existing code style
- **Update documentation:** Include documentation updates if needed

## Reporting Issues

When reporting bugs, please include:

- **Description:** Clear description of the issue
- **Steps to reproduce:** Exact steps to reproduce the problem
- **Expected behavior:** What you expected to happen
- **Actual behavior:** What actually happened
- **Environment:** R version, OS, conda environment details
- **Error messages:** Full error messages or stack traces

## Feature Requests

When suggesting enhancements:

- **Use case:** Describe why this feature would be useful
- **Proposed solution:** How you envision the feature working
- **Alternatives:** Other approaches you've considered
- **Examples:** Examples from other tools if applicable

## Code of Conduct

- Be respectful and constructive
- Welcome newcomers and help them learn
- Focus on what is best for the project
- Show empathy towards other community members

## Questions?

If you have questions about contributing, feel free to:

- Open an issue with the "question" label
- Review existing issues and discussions
- Check the README and documentation first

## License

By contributing, you agree that your contributions will be licensed under the same MIT License that covers the project.
