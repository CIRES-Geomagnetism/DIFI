# Packaging Python Projects

Upload the package to TestPyPI for developing and testing. 

## Generating distribution archives

`python3 -m pip install --upgrade build`

Now run this command from the same directory where pyproject.toml is located:
`python3 -m build`

## Uploading the distribution archives

1. Register an account at https://test.pypi.org/account/register/
2. Create one at https://test.pypi.org/manage/account/#api-tokenss, setting the “Scope” to “Entire account”.
3. Use twine to upload the distribution packages
     `python3 -m pip install --upgrade twine`
    `python3 -m twine upload --repository testpypi dist/*`
4. Once uploaded, your package should be viewable on TestPyPI. For example, https://test.pypi.org/project/DIFI/

## Installing your newly uploaded package at TestPyPI

`python3 -m pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ DIFI`

