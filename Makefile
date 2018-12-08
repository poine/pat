

test_web:
	cd docs; jekyll serve

clean:
	find . -name '*~' -exec rm {} \;
	find . -name '*.pyc' -exec rm {} \;
	rm -rf docs/_site
