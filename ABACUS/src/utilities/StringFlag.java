package utilities;

public class StringFlag implements Flag<String> {

	private final String name;
	private final String description;
	private String value;
	
	public StringFlag(String name, String value, String description) {
		this.name = name;
		this.value = value;
		this.description = description;
	}
	
	@Override
	public String getName() {
		return name;
	}

	@Override
	public String getDescription() {
		return description;
	}

	@Override
	public void setValue(String value) {
		this.value = value;
	}

	@Override
	public String getValue() {
		return value;
	}
}
