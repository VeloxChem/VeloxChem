test:
	python3 -m pytest -v --cov=veloxchem --cov-report=term --cov-report=html
serve:
	cd htmlcov && python3 -c 'import webbrowser; webbrowser.open_new_tab("http://localhost:8000")' && python3 -m http.server 
