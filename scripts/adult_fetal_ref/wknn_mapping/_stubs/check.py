def check_deps(*names):
    import importlib
    for n in names:
        try:
            importlib.import_module(n if n != "cuml" else "cuml")
        except ImportError as e:
            raise ImportError(f"Optional dependency '{n}' is required: {e}") from e
