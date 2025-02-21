import struct


def read_simh_tape(file_path):
    with open(file_path, "rb") as f:
        while True:
            # Read the record length (4 bytes, little-endian)
            length_data = f.read(4)
            if not length_data:
                break  # End of file

            record_length = struct.unpack("<I", length_data)[0]

            # Handle special cases
            if record_length == 0xFFFFFFFF:
                print("EOF marker found.")
                continue
            elif record_length == 0x00000000:
                print("EOT marker found.")
                break

            # Read the data block
            data = f.read(record_length)

            # Verify the trailing record length
            trailing_length = struct.unpack("<I", f.read(4))[0]
            if trailing_length != record_length:
                raise ValueError("Record length mismatch!")

            # Process the data (example: print it as text)
            print(f"Record length: {record_length}, Data: {data}")


# Usage
read_simh_tape("1/result.000")
