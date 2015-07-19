package utilities;

public class IntegerFlag implements Flag<Integer> {

	private final String name;
	private final String description;
	private Integer value;
	
	public IntegerFlag(String name, Integer value, String description) {
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
	public void setValue(Integer value) {
		this.value = value;
	}

	@Override
	public Integer getValue() {
		return value;
	}
}
