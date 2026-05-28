import qrcode

# Replace with your uploaded GIF URL
gif_url = "https://drive.google.com/file/d/1vhvTyE21G3JZsh_ksHYR7i69D8AwsRNX/view?usp=sharing"

# Generate QR code
qr = qrcode.QRCode(
    version=1,
    error_correction=qrcode.constants.ERROR_CORRECT_L,
    box_size=10,
    border=4,
)
qr.add_data(gif_url)
qr.make(fit=True)

# Save as an image
img = qr.make_image(fill="black", back_color="white")
img.save("qr_code.png")

print("QR code saved as qr_code.png")
