run:
	@for file in Data/*.tsp; do \
		timeout 30m python3 tp2.py $$file || true; \
	done
