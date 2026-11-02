fn main() {
    if std::env::args().nth(1).as_deref() == Some("full") { for row in sheaf_kit::full_table() { println!("{row}"); } }
    else { for row in sheaf_kit::bridge_table() { println!("{row}"); } }
}
