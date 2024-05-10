./lint.sh && \
./test.sh && \
python3 -m pip install --upgrade build && \
python3 -m pip install --upgrade twine && \
rm -rf dist && \
python3 -m build && \
git --version && \
python3 -m twine upload dist/* && \
git tag "$(python3 vaxrank/version.py)" &&  \
git push --tags

