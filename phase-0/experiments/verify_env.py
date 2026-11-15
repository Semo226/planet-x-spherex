"""Phase 0 environment verification.

Confirms that all required scientific packages import and report versions.
Run from the repo root with the project venv active:

    .venv\\Scripts\\python.exe phase-0/experiments/verify_env.py
"""

import sys

REQUIRED = {
    "numpy": "numpy",
    "pandas": "pandas",
    "matplotlib": "matplotlib",
    "astropy": "astropy",
    "astroquery": "astroquery",
    "pyvo": "pyvo",
    "requests": "requests",
}


def main() -> int:
    print(f"Python {sys.version.split()[0]} on {sys.platform}")
    failures = []
    for module, dist in REQUIRED.items():
        try:
            mod = __import__(module)
            version = getattr(mod, "__version__", "unknown")
            print(f"  OK   {module:12s} {version}")
        except Exception as exc:  # noqa: BLE001 - report and continue
            failures.append(module)
            print(f"  FAIL {module:12s} {exc}")

    if failures:
        print(f"\nFAILED imports: {', '.join(failures)}")
        return 1
    print("\nAll required packages import successfully.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
