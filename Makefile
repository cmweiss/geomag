help:
	@echo "Usage:"
	@echo "make help           -- display this help"
	@echo "make test           -- run the tests"
	@echo "make install        -- install for development"

install:
	pip install -r requirements.txt
	pip install -e .

test:
	pytest
