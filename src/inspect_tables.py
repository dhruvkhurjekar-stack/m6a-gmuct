from pathlib import Path
import pandas as pd

SUPP_PATH = Path("/Users/dhruvkhurjekar/biol3999-m6a-gmuct/biol3999-m6a-gmuct/data/Supplementary Tables.xlsx")
M6A_PATH = Path("/Users/dhruvkhurjekar/biol3999-m6a-gmuct/biol3999-m6a-gmuct/data/m6a_table_s1.xlsx")

def show_excel_info(path: Path, max_sheets: int = 20) -> None:
    print("\n" + "=" * 80)
    print(f"FILE: {path}")
    print("=" * 80)

    if not path.exists():
        print("File not found.")
        return

    xls = pd.ExcelFile(path)
    print(f"Sheets ({len(xls.sheet_names)}):")
    for s in xls.sheet_names[:max_sheets]:
        print("  -", s)

    # Read first sheet to show columns + preview
    first_sheet = xls.sheet_names[0]
    df = pd.read_excel(path, sheet_name=first_sheet)
    print("\nPreview of FIRST sheet:", first_sheet)
    print("Columns:", list(df.columns))
    print(df.head(5).to_string(index=False))

def main():
    show_excel_info(SUPP_PATH)
    show_excel_info(M6A_PATH)

if __name__ == "__main__":
    main()