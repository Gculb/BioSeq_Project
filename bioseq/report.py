import json
import os


def write_report(report_data, output_dir, report_name):
    """Write pipeline measurements as JSON and plain-text reports."""
    os.makedirs(output_dir, exist_ok=True)
    safe_name = os.path.basename(report_name)
    if not safe_name or safe_name in {".", ".."}:
        raise ValueError("report_name must be a non-empty file name.")
    safe_name = os.path.splitext(safe_name)[0]

    json_path = os.path.join(output_dir, f"{safe_name}.json")
    text_path = os.path.join(output_dir, f"{safe_name}.txt")

    with open(json_path, "w", encoding="utf-8") as file_handle:
        json.dump(report_data, file_handle, indent=2, sort_keys=True)
        file_handle.write("\n")

    with open(text_path, "w", encoding="utf-8") as file_handle:
        file_handle.write("BioSeq analysis report\n")
        file_handle.write("======================\n\n")
        file_handle.write(f"Input: {report_data['input']}\n")
        file_handle.write(f"Format: {report_data['format'].upper()}\n")
        record_count = report_data["validation"]["record_count"]
        if record_count is not None:
            file_handle.write(f"Records: {record_count}\n")
        if report_data.get("reference"):
            file_handle.write(f"Reference: {report_data['reference']}\n")
        file_handle.write("\n")
        file_handle.write("Measurements\n------------\n")
        for name, value in report_data["metrics"].items():
            label = name.replace("_", " ").capitalize()
            file_handle.write(f"{label}: {value}\n")

    return {"json": json_path, "text": text_path}
